%% Uses spm12 to smooth brain maps received with fMRIPrep
% Christian Paret, ZI-Mannheim, 2021

% Requires dataset according to BIDS standard and derivatives output from fmriPrep (tested with version 20.2.0)
% Requires dataset-table according to RaBIDS standard
%
% Unzip nifti-images before running this script.
% Download program to dataset/code directory to run

%% Read data from datasheet
data = readtable('datasheet.xlsx','ReadRowNames',true,'PreserveVariableNames',true,'NumHeaderLines',0);

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
reqtask = input('What task to work on? Enter a single task name or enter ''all'' to process all tasks available.\n');

%% Smoothing kernel
smkernel = input('Define smoothing kernel (mm).\n');

%% Initialize SPM
spm_jobman('initcfg')

%% Read confounds timeseries file and write realignment parameters to nuisance file
allsubs = dir([derivd,filesep,'sub-*']);

for sub = 1:length(allsubs)
    subd = [derivd,filesep,allsubs(sub).name];
    
    if isfolder(subd)
        
        fprintf(['\n',allsubs(sub).name,'\n'])
        
        try
            allses = dir([subd,filesep,'ses-*']);
            
            for ses = 1:length(allses)
                sesd = [subd,filesep,allses(ses).name];
                fprintf([allses(ses).name,'\n'])
                
                if isfolder([sesd,filesep,'func'])
                    niif = dir([sesd,filesep,'func',filesep,'sub-*-preproc_bold.nii']);
                    
                    for funcs = 1:length(niif)
                        pos1 = strfind(niif(funcs).name,'task');
                        pos2 = strfind(niif(funcs).name,'_space');
                        taskid = niif(funcs).name(pos1+5:pos2-1);
                        
                        if strcmp(reqtask,'all') || strcmp(reqtask,taskid)
                            pos = strfind(niif(funcs).name,'_bold.nii');
                            newnii_name = [niif(funcs).name(1:pos),'desc-s',num2str(smkernel),niif(funcs).name(pos:end)];
                            
                            if ~isfile(fullfile(derivd,allsubs(sub).name,'ses-post','func',newnii_name))
                                getscans = cellstr(spm_select('Expand',fullfile(derivd,allsubs(sub).name,'ses-post','func',niif(funcs).name)));
                                matlabbatch{1}.spm.spatial.smooth.data = getscans;

                                matlabbatch{1}.spm.spatial.smooth.fwhm = [smkernel smkernel smkernel];
                                matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                                matlabbatch{1}.spm.spatial.smooth.im = 0;
                                matlabbatch{1}.spm.spatial.smooth.prefix = 's';
                                
                                fprintf(['Start smoothing of ',niif(funcs).name,'.\n'])
                                spm_jobman('run', matlabbatch);
                                clear matlabbatch

                                fprintf('Rename to BIDS standard.\n')
                                movefile(fullfile(derivd,allsubs(sub).name,'ses-post','func',['s',niif(funcs).name]),fullfile(derivd,allsubs(sub).name,'ses-post','func',newnii_name))
                                
                                fprintf(['Saved to file ',newnii_name,'.\n'])                                
                            else
                                fprintf('Found existing nifti file for this task with same smoothing kernel. Continue with next subject/session.\n')
                            end
                            
                        end
                    end
                else
                    fprintf('No func directory found or nifti file not readable.\n');
                end
            end
        catch
            if isfolder([subd,filesep,'func'])
                niif = dir([subd,filesep,'func',filesep,'sub-*-confounds_timeseries.tsv']);
                
                for funcs = 1:length(niif)
                    pos1 = strfind(niif(funcs).name,'task');
                    pos2 = strfind(niif(funcs).name,'_desc-confounds_timeseries');
                    taskid = niif(funcs).name(pos1+5:pos2-1);
                    
                    if strcmp(reqtask,'all') || strcmp(reqtask,taskid)
                        getscans = cellstr(spm_select('Expand',fullfile(derivd,allsubs(sub).name,'ses-post','func',niif(funcs).name)));
                        matlabbatch{1}.spm.spatial.smooth.data = getscans;
                        
                        matlabbatch{1}.spm.spatial.smooth.fwhm = [smkernel smkernel smkernel];
                        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                        matlabbatch{1}.spm.spatial.smooth.im = 0;
                        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
                        
                        spm_jobman('initcfg')
                        fprintf(['Start smoothing of ',niif(funcs).name,'.\n'])
                        spm_jobman('run', matlabbatch);
                        
                        pos = strfind(niif(funcs).name,'_bold.nii');
                        newnii_name = [niif(funcs).name(1:pos),'_desc-s',num2str(smkernel),niif(funcs).name(pos:end)];
                        
                        fprintf('Rename to BIDS standard.\n')
                        movefile(fullfile(derivd,allsubs(sub).name,'ses-post','func',['s',niif(funcs).name]),fullfile(derivd,allsubs(sub).name,'ses-post','func',newnii_name))
                    end
                end
            else
                fprintf('No func directory found or nifti file not readable.\n');
            end
        end
    end
end
fprintf('End\n')