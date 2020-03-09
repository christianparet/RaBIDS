function out = firstlevel(subject,task,data_analysis_path,ses_id,addsub,condfile,overwrite)
% v0.1 release
% Estimate spm first-level model with spm12

%% Need to check whether firstlevel.m works properly because accidently overwritten by PM modulation stuff, and then deleted it to bring it byk in old shape (hopefully all fine)

%% Comment out if function in use
% clear
% clc
% data_analysis_path = 'S:\AG-Austausch\RaBIDS\example\your project directory';
% subject = 'RaBIDS01';
% ses_id = 'ses-01';
% addsub = 'y'; % add prefix 'sub-' before subject name
% task = 'faces';
% overwrite = 'y';
% condfile = 'S:\AG-Austausch\RaBIDS\example\your project directory\dataset\code\conditions_faces.xlsx';
%%

%% Define contrasts
condata = readtable(condfile,'ReadRowNames',true);

Namecol = find(strcmp(condata.Properties.VariableNames,'Name'));
% ContrastTypecol = find(strcmp(condata.Properties.VariableNames,'Contrast type')); % not yet implemented
ContrastPlus1col = find(strcmp(condata.Properties.VariableNames,'ContrastPlus1'));
ContrastMinus1col = find(strcmp(condata.Properties.VariableNames,'ContrastMinus1'));

conlines = find(contains(condata.Properties.RowNames,'Contrast'));

%% Define path
if contains(addsub,'y')
    prefix = 'sub-';
else
    prefix = '';
end

if strcmp(ses_id,'none')
    write_ses = '_';
    wd = [prefix,subject];
else
    write_ses = ['_',ses_id,'_'];
    wd = [prefix,subject,filesep,ses_id];
end

derivp = fullfile(data_analysis_path,'dataset','derivatives\RaBIDS-prepro\',wd);
sotp = fullfile(data_analysis_path,'dataset',wd,'func');
firstlevelp = fullfile(data_analysis_path,'dataset','derivatives\RaBIDS-firstlevel\',['task-',task],wd);

if ~isfolder(firstlevelp) || contains(overwrite,'y')

    %% Define model
    matlabbatch{1}.spm.stats.fmri_spec.dir = {firstlevelp}; 
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; 
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2; 
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 36; 
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 18; 

    try
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('expand',fullfile(derivp,['swa',prefix,subject,write_ses,'task-',task,'_bold.nii'])));
    catch
        out = 'Images not found. Continue with next session.\n';
        return
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {fullfile(sotp,[prefix,subject,write_ses,'task-',task,'_multicond.mat'])};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {}); 

    try
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(derivp,['nuisance_',prefix,subject,write_ses,'task-',task,'_bold.mat'])};
    catch
        out = 'Nuisance file not foun. Continue with next session.\n';
        return
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128; 
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {}); 
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; 
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1; 
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None'; 
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8; 
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''}; 
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; 

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

    % Define contrasts
    for i = 1:length(conlines)
        try
            load(matlabbatch{1}.spm.stats.fmri_spec.sess.multi{:}); % load file defining response function (i.e. stimulus onset times)
        catch
            out = 'Condition definition file not found. Continue with next session.\n';
            return
        end
        
        weights = [];
        for j = 1:length(names)
            if strcmp(condata{conlines(i),ContrastPlus1col}{:},names{j})
                weights = [weights 1];
            elseif strcmp(condata{conlines(i),ContrastMinus1col}{:},names{j})
                weights = [weights -1];
            else
                weights = [weights 0];
            end
        end
                
        if ~any(weights)
            out=['Contrast ',condata{conlines(i),Namecol}{:},' is invalid. Program stops.\n'];
            return
        end
        eval(['matlabbatch{3}.spm.stats.con.consess{',num2str(i),'}.tcon.name = ''',condata{conlines(i),Namecol}{:},''';']);
        eval(['matlabbatch{3}.spm.stats.con.consess{',num2str(i),'}.tcon.weights = ',mat2str(weights),';']);
        eval(['matlabbatch{3}.spm.stats.con.consess{',num2str(i),'}.tcon.sessrep = ''none'';']);
    end
        
    matlabbatch{3}.spm.stats.con.delete = 1;

    %% Estimate model and return
    spm_jobman('run', matlabbatch);
    out = 'Model estimation successful!\n';
    
else
    out = 'Firstlevel directory exists, continue with next session.\n';
end