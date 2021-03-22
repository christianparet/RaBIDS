%% Artifact detection and repair pipeline: Make a new directory that contains the preprocessed derivative images in 3d-nifti format
% This is necessary to run ARTREPAIR toolbox (https://cibsr.stanford.edu/tools/human-brain-project/artrepair-software.html)
% Christian Paret, ZI-Mannheim, 2021

% Requires dataset according to BIDS standard and derivatives output from fmriPrep (tested with version 20.2.0)
% Requires dataset-table according to RaBIDS standard
% FMRIPrep output images (i.e., zipped 4d-nii files) need to be unzipped first (use unzip_niftis.m to unzip derivative files imported with RaBIDS, see auxiliary files on github directory)
%
% Download program to dataset/code directory to run this script
% Make sure that ARTREPAIR is installed as spm12 toolbox, see instructions on artrepair homepage
%%
clear
clc

%% Provide information about your imaging data here
smkernel = '6'; % smoothing kernel of functional images

%% Read data from datasheet
data = readtable('datasheet.xlsx','ReadRowNames',true,'PreserveVariableNames',true,'NumHeaderLines',0);
userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));
DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};
datasetd = [data_analysis_path,filesep,'dataset'];
first_imageline = find(strcmp(data.Properties.RowNames,'first image'));
if contains(data{first_imageline,userInputcol}{:},'y')
    first_image = 1; % X=1:[first_image-1] images are skipped in create sots step and were deleted from the dicom-directory (the delete-mode will be removed in a future RaBIDS version and should not be used
    fimgs_exist = 0;
else
    first_image = data{first_imageline,minImagescol}; % X=1:[first_image-1] images are skipped in create sots step
    fimgs_exist = 1;
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

%% Define contrasts
condfile = ['conditions_',reqtask,'.xlsx'];
condata = readtable(fullfile(datasetd,'code',condfile),'ReadRowNames',true,'ReadVariableNames',true,'NumHeaderLines',0);

Namecol = find(strcmp(condata.Properties.VariableNames,'Name'));
% ContrastTypecol = find(strcmp(condata.Properties.VariableNames,'Contrast type')); % not yet implemented
ContrastPlus1col = find(strcmp(condata.Properties.VariableNames,'ContrastPlus1'));
ContrastMinus1col = find(strcmp(condata.Properties.VariableNames,'ContrastMinus1'));

conlines = find(contains(condata.Properties.RowNames,'Contrast'));

%% Go through subjects and create 3d-nifti data
allsubs = dir([derivd,filesep,'sub-EFP*']);

for subject = 1:length(allsubs)
    if allsubs(subject).isdir
        allses = dir(fullfile(derivd,allsubs(subject).name,'ses-*'));
        
        for session = 1:length(allses)
            if allses(session).isdir
                fprintf(['\nProcessing ',allsubs(subject).name,' ',allses(session).name,'.\n']);
                
                % Path to save 3d-nii data
                threedniip = fullfile(data_analysis_path,'spm_analysis','preprocessed_3dnii',allsubs(subject).name,allses(session).name,['task-',reqtask]);                

                % Search smoothed derivative data
                derivname = [allsubs(subject).name,'_',allses(session).name,'_task-',reqtask,'_*desc-preproc_desc-s',smkernel,'_bold.nii'];
                derivp = fullfile(derivd,allsubs(subject).name,allses(session).name,'func');
                derivniif = dir(fullfile(derivp,derivname));
                if length(derivniif) == 1
                    fprintf(['Found derivative file ',derivniif.name,'.\n']);
                    deriv_ok = 1;
                elseif isempty(derivniif)
                    fprintf('Derivative image not found. Skip session.\n');
                    deriv_ok = 0;
                else
                    fprintf('Too many matching derivative files. Skip session.\n');
                    deriv_ok = 0;
                end
                
                if deriv_ok
                    % Check whether 3d-nii data already exist
                    pos = strfind(derivniif.name,'_bold');
                    threedniif = dir(fullfile(threedniip,[derivniif.name(1:pos),'bold_0001.nii']));
                    if length(threedniif) == 1
                        fprintf('Found at least one 3d-nifti file in target directory. Skip session.\n');
                        threednii_ok = 0;
                    else
                        threednii_ok = 1;
                        if ~exist(threedniip,'dir') % if directory not exists yet: create output dir
                            mkdir(threedniip)
                        end
                    end
                else
                    threednii_ok = 0;
                end
                
                if deriv_ok && threednii_ok
                    % Code comes from Chris Roden, https://www.nitrc.org/forum/forum.php?thread_id=9969&forum_id=4703
                    [pth,nam,ext] = spm_fileparts(fullfile(derivp,derivniif.name));
                    image = fullfile(pth,[nam,ext]); %'img.nii,1' -> 'img.nii'
                    hdr = spm_vol(image);  
                    fprintf('Reading volumes...\n')
                    img = spm_read_vols(hdr);
                    nvol = numel(hdr);
                    hdr = hdr(1);
                    fprintf('Writing volumes to output directory...\n')
                    for vol=1: nvol
                        hdr.fname = fullfile(pth, [nam, '_', num2str(vol,'%04d'), ext]);
                        spm_write_vol(hdr,img(:, :, :, vol));
                        movefile(fullfile(pth, [nam, '_', num2str(vol,'%04d'), ext]),fullfile(threedniip, [nam, '_', num2str(vol,'%04d'), ext]))
                    end
                end
            end
        end
    end
end

fprintf('End.\n')