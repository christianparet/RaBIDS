%% Estimate spm firstlevel with spm12
% Christian Paret, ZI-Mannheim, 2021

% Requires dataset according to BIDS standard and derivatives output from fmriPrep (tested with version 20.2.0)
% Requires dataset-table according to RaBIDS standard
%
% You should create nuisance-regressor files first. The program "create_nuiscnace_regfile_for_spm" (see RaBIDS/auxiliary tools) writes nuisance files to the derivatives folder that are automatically used by this script.
% You should also provide a (unzipped) smoothed, preprocessed image in the derivatives directory. Use the "smooth_images.m" program (see RaBIDS\auxiliary tools); it will automatically provide the images in the needed naming format.
%
% Download program to dataset/code directory to run
%%
clear
clc

%% Provide information about your imaging data here:
smkernel = '6'; % smoothing kernel of functional images
TR = 2; % Repetition Time
slices = 36; % number of slices per volume

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

%% Go through subjects and calculate firstlevel
allsubs = dir([derivd,filesep,'sub-EFP*']);

for subject = 3%1:length(allsubs)
    if allsubs(subject).isdir
        allses = dir(fullfile(derivd,allsubs(subject).name,'ses-*'));
        
        for session = 1:length(allses)
            if allses(session).isdir
                fprintf(['\nProcessing ',allsubs(subject).name,' ',allses(session).name,'.\n']);
                
                %% Define paths and do couple of checks before creating firstlevel model
                % Is there already an SPM.mat?
                firstleveld = [allsubs(subject).name,'_',allses(session).name,'_task-',reqtask];
                firstlevelp = fullfile(data_analysis_path,'spm_analysis','firstlevel',allses(session).name,['task-',reqtask],'taskrelated_activity',firstleveld);
                if ~isfile([firstlevelp,filesep,'SPM.mat'])
                    spm_ok = 1;
                else
                    fprintf('Model appears to be estimated already, as an SPM.mat exists in the target directory. Skip session.\n');
                    spm_ok = 0;
                end
                
                % Does the nuisance regressor file exist?
                multiregf = [allsubs(subject).name,'_',allses(session).name,'_task-',reqtask,'_desc-confounds_timeseries_desc-realignment_regressors.mat'];
                multiregp = fullfile(derivd,allsubs(subject).name,allses(session).name,'func',multiregf);
                if isfile(multiregp)
                    multireg_ok = 1;
                else
                    fprintf('Nuisance file not found. Skip session.\n');
                    multireg_ok = 0;
                end
                
                % Is there a Stimulus Onset Times file defining the stimulus protocol?
                multicondf = [allsubs(subject).name,'_',allses(session).name,'_task-',reqtask,'_multicond.mat'];
                multicondp = fullfile(datasetd,allsubs(subject).name,allses(session).name,'func',multicondf);
                if isfile(multicondp)
                    multicond_ok = 1;
                else
                    fprintf('SOTs file not found. Skip session.\n');
                    multicond_ok = 0;
                end
                
                % Does a smoothed derivative image exist?
                derivname = [allsubs(subject).name,'_',allses(session).name,'_task-',reqtask,'_*desc-preproc_desc-s',smkernel,'_bold.nii'];
                derivp = fullfile(derivd,allsubs(subject).name,allses(session).name,'func',derivname);
                derivnii = dir(derivp);
                if length(derivnii) == 1
                    fprintf(['Found derivative file ',derivnii.name,'.\n']);
                    deriv_ok = 1;
                elseif isempty(derivnii)
                    fprintf('Derivative image not found. Skip session.\n');
                    deriv_ok = 0;
                else
                    fprintf('Too many matching derivative files. Skip session.\n');
                    deriv_ok = 0;
                end
                
                if spm_ok && multireg_ok && multicond_ok && deriv_ok
                    %% Make and execute batch
                    % Define spm12 model
                    matlabbatch{1}.spm.stats.fmri_spec.dir = {firstlevelp};
                    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
                    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = slices;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = round(slices/2); % Input reference slice to which images were aligned to during slice-timing step. fMRIPrep aligns to middle slice if corresponding metadata information is available.
                    
                    getscans = cellstr(spm_select('Expand',fullfile(derivd,allsubs(subject).name,allses(session).name,'func',derivnii.name)));
                    
                    if fimgs_exist
                        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = getscans(first_image:end);
                    else
                        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = getscans;
                    end
                    
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {multicondp};
                    
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {multiregp};
                    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
                    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
                    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
                    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
                    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.3;
                    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
                    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
                    
                    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
                    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
                    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
                    
                    multiconddata = load(matlabbatch{1}.spm.stats.fmri_spec.sess.multi{:});
                    
                    for i = 1:length(conlines)
                        weights = [];
                        for j = 1:length(multiconddata.names)
                            if strcmp(condata{conlines(i),ContrastPlus1col}{:},multiconddata.names{j})
                                weights = [weights 1];
                            elseif strcmp(condata{conlines(i),ContrastMinus1col}{:},multiconddata.names{j})
                                weights = [weights -1];
                            else
                                weights = [weights 0];
                            end
                        end
                        
                        if ~any(weights)
                            fprintf(['Contrast ',condata{conlines(i),Namecol}{:},' is invalid.\n']);
                            return
                        end
                        eval(['matlabbatch{3}.spm.stats.con.consess{',num2str(i),'}.tcon.name = ''',condata{conlines(i),Namecol}{:},''';']);
                        eval(['matlabbatch{3}.spm.stats.con.consess{',num2str(i),'}.tcon.weights = ',mat2str(weights),';']);
                        eval(['matlabbatch{3}.spm.stats.con.consess{',num2str(i),'}.tcon.sessrep = ''none'';']);
                    end
                    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
                    matlabbatch{3}.spm.stats.con.delete = 1;
                    
                    %% Estimate model
                    spm_jobman('run', matlabbatch);
                    fprintf('Model estimation successful!\n');
                    clear matlabbatch
                end
            end
        end
    end
end

fprintf('End.\n')