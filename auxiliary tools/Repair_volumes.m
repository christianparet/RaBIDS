%% Artifact detection and repair pipeline: Make a new directory with the repaired, preprocessed derivative images in 3d-nifti format
% This is a wrapper around art_global which is part of the ArtRepair toolbox for spm, which can be downloaded here: https://cibsr.stanford.edu/tools/human-brain-project/artrepair-software.html
% Has been tested with spm12 and ArtRepair version 5b
% Make sure that ArtRepair is installed as spm12 toolbox, see instructions on artrepair homepage
%
% Run this script to have data ready to run art_redo
% Christian Paret, ZI-Mannheim, 2021

% Requires dataset according to BIDS standard and derivatives output from fmriPrep (tested with version 20.2.0)
% Requires dataset-table according to RaBIDS standard
% FMRIPrep output images (i.e., zipped 4d-nii files) need to be unzipped first (use unzip_niftis.m to unzip derivative files imported with RaBIDS, see auxiliary files on github directory)
%
% Note: We use the datasheet entries to ObjectType first image to select
% images for repair. Only images are imported that will go into the SPM GLM
% model. This is necessary because ArtRepair v5b reports repaired and
% deweighted images relative to the images submitted to the art_global.m
% function. Thus, if the user was going to use a subsample of repaired
% images (for art_redo.m) after running art_global.m, the text-files including the repaired and deweighted
% images would be shifted by the number of the discarded files.
%
% Download program to dataset/code directory to run this script

%% Change log
% 2021/08/03: to identify tsv-file including regressors, search word "confounds" was changed to wild card for  downwards compatibility with fMRIPrep versions before 20.2.0
%%

clear
clc

fprintf('Repair volumes using ArtRepair.\n\n')

%% Provide information about your imaging data here
smkernel = '6'; % smoothing kernel of functional images; program looks for derivatives with this smoothing kernel
use_trimmed = 1; % default = 1. program searches for a trimsession file in the derivatives directory with format <subject ID>_<session ID>_task-<task ID>_trimsession.mat. The file must have two entries pointing a start volume and a stop volume for trimming the session data.
mv_threshold = 1.0 ; % recommendation from ArtRepair toolbox (see art_global header information): good data: 0.3, moderately noisy data: 0.5, severely noisy data: 1.0

%% Set path to ArtRepair toolbox
addpath('C:\Program Files\spm12\toolbox\ArtRepair')

%% Read data from datasheet
data = readtable('datasheet.xlsx','ReadRowNames',true,'PreserveVariableNames',true,'NumHeaderLines',0);
userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));
DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};
datasetd = [data_analysis_path,filesep,'dataset'];
first_imageline = find(strcmp(data.Properties.RowNames,'first image'));
if contains(data{first_imageline,userInputcol}{:},'y')
    first_vol = 1; % X=1:[first_image-1] images were deleted from the dicom-directory (the delete-mode will be removed in a future RaBIDS version and should not be used)
    fprintf('Note: Program assumes that initial images have been deleted.\nFor more information see RaBIDS manual, ObjectType first image.\n\n')
else
    first_vol = data{first_imageline,minImagescol}; % X=1:[first_image-1] images are skipped below
    fprintf('Note: Program assumes that initial images have not been deleted.\nInitial volumes will be discarded and nuisance regressors will be cut accordingly.\nFor more information see RaBIDS manual, ObjectType first image.\n\n')
end

%% Identify derivative directory to work on
derivd = dir([data_analysis_path,filesep,'derivatives',filesep,'fmriprep*']);

fprintf(['Found ',num2str(length(derivd)),' derivative directories:\n'])

for i=1:length(derivd)
    fprintf([num2str(i),') ',derivd(i).name,'\n'])
end

try
    [c,~] = size(derivd);
catch
    fprintf('Cannot read derivatives directory.\n')
    return
end

if c > 1
    selectd = input('Enter number of derivative directory to work on and press enter.\n');
else
    selectd = 1;
end
derivd = [data_analysis_path,filesep,'derivatives',filesep,derivd(selectd).name];

% in case of nested structure: account for 2 layers of fmriprep directories
try
    childd = dir([derivd,filesep,'fmriprep*']);
    if length(childd)>1
        fprintf('Too many subdirectories.\n')
        return
    else
        derivd = [derivd,filesep,childd.name];
    end
end

try
    childd = dir([derivd,filesep,'fmriprep*']);
    if length(childd)>1
        fprintf('Too many subdirectories.\n')
        return
    else
        derivd = [derivd,filesep,childd.name];
    end
end

%% Task to work on
reqtask = input('What task to work on? Enter a single task name.\n');

%% Go through subjects
allsubs = dir([derivd,filesep,'sub-EFP*']);

for subject = 1:length(allsubs)
    if allsubs(subject).isdir
        allses = dir(fullfile(derivd,allsubs(subject).name,'ses-*'));
        
        for session = 1:length(allses)
            if allses(session).isdir
                fprintf(['\nProcessing ',allsubs(subject).name,' ',allses(session).name,'.\n']);
                
                %% Check for smoothed derivative data
                derivname = [allsubs(subject).name,'_',allses(session).name,'_task-',reqtask,'_*desc-preproc_desc-s',smkernel,'_bold.nii'];
                derivp = fullfile(derivd,allsubs(subject).name,allses(session).name,'func');
                derivniif = dir(fullfile(derivp,derivname));
                if length(derivniif) == 1
                    fprintf(['    Found derivative file ',derivniif.name,'.\n']);
                    deriv_ok = 1;
                elseif isempty(derivniif)
                    fprintf('    Derivative image not found. Skip session.\n');
                    deriv_ok = 0;
                else
                    fprintf('    Too many matching derivative files. Skip session.\n');
                    deriv_ok = 0;
                end
                
                %% Trim session?
                trim = false; % set flag false as standard. Do now couple of checks to give detailed user feedback
                start_vol = first_vol;
                trimf = [allsubs(subject).name,'_',allses(session).name,'_task-',reqtask,'_trimsession.mat'];
                threedniip = fullfile(data_analysis_path,'spm_analysis','repaired_preprocessed_3dnii',allses(session).name,allsubs(subject).name,['task-',reqtask],'complete_session_data'); % Path to save repaired 3d-nii data
                if use_trimmed && isfile(fullfile(derivp,trimf))
                    fprintf('    Trimsession file found, importing trimmed session data.\n')
                    trim = true;
                    load([derivp,filesep,trimf]) % 1-by-2 vector with information of first/last scan to use; see script "find_spikes_and_trim_session.m"
                    if start_vol < trimses(1)
                        start_vol = trimses(1);
                    end
                    threedniip = fullfile(data_analysis_path,'spm_analysis','repaired_preprocessed_3dnii',allses(session).name,allsubs(subject).name,['task-',reqtask],'trimmed_session_data'); % Path to save 3d-nii data
                elseif use_trimmed && ~isfile(fullfile(derivp,trimf))
                    fprintf('    No trimsession file found, importing full session data.\n')
                else
                    fprintf('   Importing full session data.\n')
                end
                
                %% Check for confounds time series file and read nuisance regressors into matrix R
                conf_timeseriesf = dir(fullfile(derivp,'sub-*-confounds_*.tsv')); 
                
                if ~isempty(conf_timeseriesf)
                    fprintf('    Confound timeseries file found.\n');
                    conf_ok = 1;
                    for funcs = 1:length(conf_timeseriesf)
                        pos1 = strfind(conf_timeseriesf(funcs).name,'task');
                        pos2 = strfind(conf_timeseriesf(funcs).name,'_desc-confounds');
                        taskid = conf_timeseriesf(funcs).name(pos1+5:pos2-1);
                        
                        if strcmp(reqtask,taskid)
                            conf_table = readtable(fullfile(derivp,conf_timeseriesf(funcs).name),'FileType','text','ReadRowNames',false,'PreserveVariableNames',true,'NumHeaderLines',0);
                            
                            % realignment regressors
                            R = [conf_table.trans_x, conf_table.trans_y, conf_table.trans_z, conf_table.rot_x, conf_table.rot_y, conf_table.rot_z];
                            
                            % discard initial volumes and trim session data if needed
                            if trim
                                R = R(start_vol:trimses(2),:);
                            else
                                R = R(start_vol:end,:);
                            end
                            
                        end
                    end
                else
                    fprintf('    No confound-timeseries file found. Skip session.\n')
                    conf_ok = 0;
                end
                               
                %% Process session data if requirements are met
                if deriv_ok && conf_ok
                    
                    %% Check whether repaired 3d-nii data exists (note: ArtRepair v5b is incompatible with 4d-nifti data)
                    pos = strfind(derivniif.name,'_bold');
                    threedniif = dir(fullfile(threedniip,['v',derivniif.name(1:pos),'bold_',num2str(start_vol,'%.4d'),'.nii']));
                    if length(threedniif) == 1
                        fprintf('    Found at least one repaired volume in target directory. Skip Session.\n');
                    else
                        threednii_ok = 0;
                        if ~exist(threedniip,'dir') % if directory not exists yet: create output dir
                            mkdir(threedniip)
                        end
                        fprintf('    Convert 4d-nii to 3d-nii.\n')
                        
                        %% Create 3d-nifti data, because ArtRepair v5b is incompatible with 4d-nifti
                        % Code comes from Chris Roden, https://www.nitrc.org/forum/forum.php?thread_id=9969&forum_id=4703
                        [pth,nam,ext] = spm_fileparts(fullfile(derivp,derivniif.name));
                        image = fullfile(pth,[nam,ext]); %'img.nii,1' -> 'img.nii'
                        hdr = spm_vol(image);
                        fprintf('    Reading 4d-nii volumes...\n')
                        img = spm_read_vols(hdr);
                        nvol = numel(hdr); 
                        
                        % in case the session is to be trimmed we will only import images as set in trimsession file
                        if trim
                            stop_vol = trimses(2);
                        else
                            stop_vol = nvol;
                        end
                        
                        hdr = hdr(1);
                        fprintf('   Writing 3d-nii volumes to output directory...\n')
                        for vol = start_vol : stop_vol % only import images that go into the glm
                            hdr.fname = fullfile(pth, [nam, '_', num2str(vol,'%04d'), ext]);
                            spm_write_vol(hdr,img(:, :, :, vol));
                            movefile(fullfile(pth, [nam, '_', num2str(vol,'%04d'), ext]),fullfile(threedniip, [nam, '_', num2str(vol,'%04d'), ext]))
                        end
                        
                        %% Bad volumes: detect and repair
                        
                        % Load volumes in character variable
                        getscans = spm_select('FPList',threedniip,['^',derivniif.name(1:pos),'bold_.*\.nii$']);
                        
                        % Make a copy of nuisance regressors that is readable with ArtRepair v5b
                        writematrix(R,[threedniip,filesep,'realignment_regressors'],'Delimiter','tab');
                        
                        % Repair scans using ArtRepair toolbox. Note:
                        % art_global_mvthreshold is an adapted version of
                        % art_global.m v 2.6, expecting mv_threshold as an
                        % input
                        art_global_mvthreshold(getscans, [threedniip,filesep,'realignment_regressors.txt'], 4, 2, mv_threshold) % second last input is RepairType: only repaired scans should be deweighted. This option should be used if motion regressors are used for nuisance regression in the firstlevel GLM estimation.
                        
                        % Remove the original, unrepaired 3d-niftis as they are no longer needed
                        originscans = dir(fullfile(threedniip, [nam, '_*']));
                        fprintf('    Removing original 3d-nii volumes...\n')
                        for i=1:length(originscans)
                            delete(fullfile(threedniip,originscans(i).name));
                        end
                        
                    end
                    
                end
            end
        end
    end
end

fprintf('End.\n')