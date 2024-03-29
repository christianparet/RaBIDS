%% Use spm12 to smooth brain maps received with fMRIPrep
% Christian Paret, ZI-Mannheim, 2021

% Requires dataset according to BIDS standard and derivatives output from fmriPrep (tested with version 20.2.0)
% Requires dataset-table according to RaBIDS standard
%
% Unzip nifti-images before running this script.
% Download program to dataset/code directory to run

%% Changelog

% 2021/05/20 fixed bug to handle datasets without session-subdirectories
% (i.e., with subject directories having func directories on next lower
% hierarchy.

clear

%% Read data from datasheet
data = readtable('datasheet.xlsx','ReadRowNames',true,'PreserveVariableNames',true,'NumHeaderLines',0);
AddPathline = find(contains(data.Properties.RowNames,'add path'));
userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));
for i = 1:length(AddPathline)
    addpath(data{AddPathline(i),userInputcol}{:});
end
userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));
DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};

%% Identify derivative directory to work on
derivd = dir([data_analysis_path,filesep,'derivatives',filesep,'fmriprep*']);

fprintf(['Found ',num2str(length(derivd)),' derivative directories:\n'])

for i=1:length(derivd)
    fprintf([num2str(i),') ',derivd(i).name,'\n'])
end

selectd = input('Enter number of derivative directory to unzip preporcessed niftis and press enter.\n');
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
reqtask = input('What task to work on? Enter ''TaskID'' to process a single task or/nenter ''all'' to process all functional scans available.\n');

%% Smoothing kernel
smkernel = input('Define smoothing kernel (mm).\n');

%% Initialize SPM
spm_jobman('initcfg')

%% Smooth images
allsubs = dir([derivd,filesep,'sub-*']);

for sub = 1:length(allsubs)
    subd = [derivd,filesep,allsubs(sub).name];
    
    if isfolder(subd)
        
        fprintf([allsubs(sub).name,'\n'])
        
        allses = dir([subd,filesep,'ses-*']);        
        nses = length(allses);
        if nses < 1 % Account for datasets that do not have ses subdirectory
            nses = 1;
            sesid = false;
        else
            sesid = true;
        end
        
        for ses = 1:nses
            if sesid
                sesdir = [subd,filesep,allses(ses).name];
                fprintf([allses(ses).name,'\n'])
            else
                sesdir = subd;
            end
            
            if isfolder([sesdir,filesep,'func'])
                niif = dir(fullfile(sesdir,'func','sub-*-preproc_bold.nii'));
                
                for funcs = 1:length(niif)
                    pos1 = strfind(niif(funcs).name,'task');
                    pos2 = strfind(niif(funcs).name,'_space');
                    taskid = niif(funcs).name(pos1+5:pos2-1);
                    
                    if strcmp(reqtask,'all') || strcmp(reqtask,taskid)
                        pos1 = strfind(niif(funcs).name,'_bold.nii');
                        newnii_name = [niif(funcs).name(1:pos1),'desc-s',num2str(smkernel),niif(funcs).name(pos1:end)];
                        
                        if ~isfile(fullfile(sesdir,'func',newnii_name))
                            getscans = cellstr(spm_select('Expand',fullfile(sesdir,'func',niif(funcs).name)));
                            matlabbatch{1}.spm.spatial.smooth.data = getscans;
                            
                            matlabbatch{1}.spm.spatial.smooth.fwhm = [smkernel smkernel smkernel];
                            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                            matlabbatch{1}.spm.spatial.smooth.im = 0;
                            matlabbatch{1}.spm.spatial.smooth.prefix = 's';
                            
                            fprintf(['Start smoothing of ',niif(funcs).name,'.\n'])
                            spm_jobman('run', matlabbatch);
                            clear matlabbatch
                            
                            fprintf('Rename to BIDS standard.\n')
                            movefile(fullfile(sesdir,'func',['s',niif(funcs).name]),fullfile(sesdir,'func',newnii_name))
                            
                            % add json file for smoothed image
                            pos2 = strfind(niif(funcs).name,'.nii');
                            jsonf = [niif(funcs).name(1:pos2),'json'];
                            jsondata = jsondecode(fileread(fullfile(sesdir,'func',jsonf)));
                            jsondata.SmoothingKernel = smkernel;
                            pos3 = strfind(newnii_name,'.nii');
                            newjsonf = [newnii_name(1:pos3),'json'];
                            saveJSONfile(jsondata,fullfile(sesdir,'func',newjsonf)) % Lior Kirsch (2020). Structure to JSON (https://www.mathworks.com/matlabcentral/fileexchange/50965-structure-to-json), MATLAB Central File Exchange. Retrieved February 27, 2020.
                            
                            fprintf(['Saved to file ',newnii_name,'.\n'])
                        else
                            fprintf('Found existing nifti file for this task with same smoothing kernel. Continue with next subject/session.\n')
                        end
                        
                    end
                end
            else
                fprintf('No func directory found.\n');
            end
        end
        
    end
end
fprintf('End\n')