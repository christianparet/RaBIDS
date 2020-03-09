function out = prepro(subject,task,data_analysis_path,spm_path,ses_id,addsub,overwrite)
% v0.1 release

%% Automatic preprocessing of functional fMRI data with spm12

%% Comment out if function in use
% clear
% subject= 'RaBIDS01';
% task = 'scenes'; % name task for identification of functional runs in dataset
% data_analysis_path = 'S:\AG-Austausch\RaBIDS\example\your project directory';
% spm_path = 'C:\Program Files\spm12'; % path to spm12
% ses_id = 'ses-01';
% addsub = 'yes';
% overwrite = 'y';
%%

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

% Check for existing file
derivp = fullfile(data_analysis_path,'dataset','derivatives\RaBIDS-prepro\',wd);
derivf = ['swa',prefix,subject,write_ses,'task-',task,'_bold.nii'];

if ~isfile(fullfile(derivp,derivf)) || contains(overwrite,'y')
    
    try
        matlabbatch{1}.spm.temporal.st.scans = {cellstr(spm_select('expand',fullfile(data_analysis_path,'dataset',wd,'func',[prefix,subject,write_ses,'task-',task,'_bold.nii']))) };
    catch
        out = 'Functional data not found. Continue with next session.\n';
        return
    end

    matlabbatch{1}.spm.temporal.st.nslices = 36;
    matlabbatch{1}.spm.temporal.st.tr = 2;
    matlabbatch{1}.spm.temporal.st.ta = 1.94444;
    matlabbatch{1}.spm.temporal.st.so = [36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1];
    matlabbatch{1}.spm.temporal.st.refslice = 18;
    matlabbatch{1}.spm.temporal.st.prefix = 'a';
    matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [0 1];
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    matlabbatch{3}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));

    try
        matlabbatch{3}.spm.spatial.coreg.estimate.source = cellstr(spm_select('expand',fullfile(data_analysis_path,'dataset',wd,'anat',[prefix,subject,write_ses,'T1w.nii'])));
    catch
        out = 'Anatomical data not found. Continue with next session.\n';
        return
    end

    matlabbatch{3}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    matlabbatch{4}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{4}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{4}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{4}.spm.spatial.preproc.channel.write = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(1).tpm = {[spm_path,'\tpm\TPM.nii,1']};
    matlabbatch{4}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{4}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(2).tpm = {[spm_path,'\tpm\TPM.nii,2']};
    matlabbatch{4}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{4}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(3).tpm = {[spm_path,'\tpm\TPM.nii,3']};
    matlabbatch{4}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{4}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(4).tpm = {[spm_path,'\tpm\TPM.nii,4']};
    matlabbatch{4}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{4}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(5).tpm = {[spm_path,'\tpm\TPM.nii,5']};
    matlabbatch{4}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{4}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(6).tpm = {[spm_path,'\tpm\TPM.nii,6']};
    matlabbatch{4}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{4}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{4}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{4}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{4}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{4}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{4}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{4}.spm.spatial.preproc.warp.write = [0 1];
    matlabbatch{5}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Realign: Estimate: Realigned Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','cfiles'));
    matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                              78 76 85];
    matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{5}.spm.spatial.normalise.write.woptions.prefix = 'w';
    matlabbatch{6}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{6}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{6}.spm.spatial.smooth.dtype = 0;
    matlabbatch{6}.spm.spatial.smooth.im = 0;
    matlabbatch{6}.spm.spatial.smooth.prefix = 's';

    spm_jobman('initcfg')

    preprobatch=matlabbatch;
    clear('matlabbatch')
    spm_jobman('run', preprobatch);

    % delete intermediate files created by spm that are not further needed
    funcp = fullfile(data_analysis_path,'dataset',wd,'func');
    anatp = fullfile(data_analysis_path,'dataset',wd,'anat');
        
    fprintf('\nDelete intermediate spm files from raw data directory\n')

    fi = fullfile(funcp,['a',prefix,subject,write_ses,'task-',task,'_bold.nii']);
    delete(fi)

    fi= fullfile(funcp,['a',prefix,subject,write_ses,'task-',task,'_bold.mat']);
    delete(fi)

    fi=fullfile(funcp,['meana',prefix,subject,write_ses,'task-',task,'_bold.nii']);
    delete(fi)

    fi=fullfile(funcp,['wa',prefix,subject,write_ses,'task-',task,'_bold.nii']);
    delete(fi)
    
    for tissuemap = 1:5
        fi=fullfile(anatp,['c',num2str(tissuemap),prefix,subject,write_ses,'T1w.nii']);
        delete(fi)
    end
    
    fi=fullfile(anatp,[prefix,subject,write_ses,'T1w_seg8.mat']);
    delete(fi)
    
    fi=fullfile(anatp,['y_',prefix,subject,write_ses,'T1w.nii']);
    delete(fi)

    % copy preprocessed files to derivative directory
    fprintf('\nCopy derivative data\n')
    try
        mkdir(derivp)
    catch 
        fprintf('Copy into existing directory\n')
    end

    [fp,fn,fe] = fileparts(fullfile(funcp,derivf));
    try movefile(fullfile(fp,[fn,fe]),fullfile(derivp,[fn,fe]))
    catch
        fprintf(['File ',[fn,fe],' not found in raw data folder\n'])
    end

    [fp,fn,fe] = fileparts(fullfile(funcp,['rp_a',prefix,subject,write_ses,'task-',task,'_bold.txt']));
    try movefile(fullfile(fp,[fn,fe]),fullfile(derivp,[fn,fe]))
        catch
        fprintf(['File ',[fn,fe],' not found in raw data folder\n'])
    end

    out = 'Finished preprocessing\n';
    
else
    out = 'Derivative file exists, continue with next step.\n';
end