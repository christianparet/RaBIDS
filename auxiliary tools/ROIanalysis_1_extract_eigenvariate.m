%% Extract ROI data using SPM's VOI utility
% Christian Paret, ZI Mannheim, 2021-2022

% Script can be used on firstlevel files produced with SPManalysis_2_firstlevel.m
% Uses VOI uitilities of SPM12

%% Cave:
% During extraction of eigenvariate (i.e., the BOLD signal time course), the time course data is adjusted for the
% contrast of interest (default: contrast image #1). Usually, this is an effects of interest (eoi) F-contrast. Note that all effects
% (i.e., regressors) not included in the contrast of interest are regressed out from the time course data.
%%

clc
clear

%% Define directories and initiate variables
fail_counter = 1;
fprintf('Select directory containing subject directories with SPM.mat files.\nSomething like %s.','\mytrainingdata\your project directory\spm_analysis\firstlevel\ses-01\task-faces\taskrelated_activity\original')
firstleveldir = uigetdir(pwd,'Select firstlevel directory that includes the subject directories. See command window for help.'); % directory containing child directories with each child containing an estimated firstlevel SPM.mat
% firstleveldir = 'Y:\Projects\EFPTest\Data_analysis\spm_analysis\firstlevel\ses-pre\task-efptest\taskrelated_activity\repaired';
subdirs = dir([firstleveldir,filesep,'sub*']);
eoi = 1; % number of contrast to adjust for
%%
if ~isempty(subdirs)
    
    % To select all masks within directory:
    maskdir = uigetdir(pwd,'Select the directory including your mask files (*.nii). All contained masks will be processed.');
    maskfs = dir([maskdir,'/*.nii']);
    [maskn,~] = size(maskfs);
    
    %% Iterate through masks
    for maskcount = 1:maskn
        [~,maskf,fe] = fileparts(fullfile(maskdir,maskfs(maskcount).name));  % this is a bit odd but necessary for compatibility with the two mask-selection options above. Use if masks-directory option is used; comment out if single-mask selection is used
        fprintf(['\n------------ Processing mask: ',maskf,' ------------.\n'])
        
        %% Iterate through subjects
        for sub = 1:length(subdirs)
            if subdirs(sub).isdir
                fprintf(['\nProcessing ',subdirs(sub).name,'.\n'])
                dum1 = strfind(subdirs(sub).name,'_ses-');
                if isempty(dum1)
                    dum1 = strfind(subdirs(sub).name,'_task-');
                    if isempty(dum1)
                        fprintf('An error occured. Directory may not be consistent with BIDS.\n')
                        return
                    end
                end
                
                spm_name = fullfile(firstleveldir,subdirs(sub).name,'SPM.mat');
                
                if exist(spm_name,'file')
                    
                    spmf = load(spm_name);
                    nr_sess = length(spmf.SPM.Sess);
                    
                    for sess = 1:nr_sess
                        
                        clear matlabbatch
                        matlabbatch{1}.spm.util.voi.spmmat = {spm_name}; % path to SPM.mat
                        matlabbatch{1}.spm.util.voi.adjust = eoi; % Adjust for F Contrast of Interest; it is necessary that this number matches the order of contrasts as they were defined in the firstlevel SPM-batch.
                        matlabbatch{1}.spm.util.voi.session = sess;
                        matlabbatch{1}.spm.util.voi.name = ['mask-',maskf,'_sess']; % name to save; spm appends '_SESSION-NUMBER' to file
                        matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {fullfile(maskdir,[maskf,fe,',1'])}; % mask images defined above
                        matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5; % threshold for mask (for binary masks anything >0 % <1 is fine)
                        matlabbatch{1}.spm.util.voi.expression = 'i1';
                        
                        try
                            spm_jobman('run', matlabbatch)
                        catch
                            failure{fail_counter,1}=subdirs(sub).name; failure{fail_counter,sess}=['sess-',int2str(sess),' failed'];
                            fail_counter=fail_counter+1;
                        end
                    end
                    
                else
                    fprintf('   No SPM.mat found for this subject.\n')
                end
                
            end
            
        end
    end
end