%% Use spm12 to smooth brain maps received with fMRIPrep
% Christian Paret, ZI-Mannheim, 2021-2026

% Requires dataset according to BIDS standard and derivatives output from fmriPrep (tested with version 20.2.0)
% Requires dataset-table according to RaBIDS standard
%
% Unzip nifti-images before running this script.
% Download program to dataset/code directory to run

%% Changelog

% 2026/01/05 derivatives selection using 1-click solution

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
derivd = uigetdir(pwd,'Select subject directory with derivatives (''//derivatives/sub-XXX/sub-XXX''.');

%% Task to work on
reqtask = input('What task to work on? Enter ''TaskID'' to process a single task or/nenter ''all'' to process all functional scans available.\n');

%% Smoothing kernel
smkernel = input('Define smoothing kernel (mm).\n');

%% Initialize SPM
fprintf('Initializing SPM...\n')
% spm_jobman('initcfg')

%% Smooth images
allses = dir([derivd,filesep,'ses-*']);
nses = length(allses);
if nses < 1 % Account for datasets that do not have ses subdirectory
    nses = 1;
    sesid = false;
else
    sesid = true;
end

for ses = 1:nses
    if sesid
        sesdir = [derivd,filesep,allses(ses).name];
        fprintf([allses(ses).name,'\n'])
    else
        sesdir = derivd;
    end

    if isfolder([sesdir,filesep,'func'])
        niif = dir(fullfile(sesdir,'func','sub-*-preproc_bold.nii'));

        for funcs = 1:length(niif)
            pos1 = strfind(niif(funcs).name,'task');
            pos2 = strfind(niif(funcs).name,'_run'); % if there is no 'run' defined in task-name, '_space' should be used instead
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

                    % Create new json sidecar file
                    pos1 = strfind(niif(funcs).name,'.nii');
                    json_info = jsondecode(fileread(fullfile(sesdir,'func',[niif(funcs).name(1:pos1),'json'])));
                    pos1 = strfind(sesdir,'sub-');
                    if pos1>1
                        pos1 = pos1(2);
                    end
                    json_info.Sources{end+1} = ['bids::',sesdir(pos1:end),'/',niif(funcs).name];
                    dum1 = strfind(json_info.Sources{end},'\');
                    json_info.Sources{end}(dum1) = '/';
                    json_info.SmoothingKernelFWHM = smkernel;
                    json_info.SmoothingKernelUnit = 'mm';
                    json_info.SmoothingSoftware = spm('version');
                    json_new = jsonencode(json_info,PrettyPrint=true);
                    pos1 = strfind(newnii_name,'.nii');
                    fid = fopen(fullfile(sesdir,'func',[newnii_name(1:pos1),'json']), 'w');
                    fprintf(fid, '%s', json_new);
                    fclose(fid);

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