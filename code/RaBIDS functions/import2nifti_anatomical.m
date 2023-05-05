function out = import2nifti_anatomical(HowExpectDicoms,dicomdir,subject,suff,ses_id,study_identifier,data_analysis_path,series,series_n,addsub,overwrite)
% Most recently updated for v0.2.2 release

% Uses dicm2nii toolbox to import dicom-images to nifti-format
% Folder structure in BIDS format is produced, for more information see: Gorgolewski, K. J. et al., Sci. Data 3:160044, doi: 10.1038/sdata.2016.44 (2016)

% Christian Paret, Central Institute of Mental Health Mannheim, 2018-2020

% Change log from v0.1:
% - implemented json-output check
% - implemented Warning-reporting
% - account for nii.gz

%% Comment out if function in use
% HowExpectDicoms = 'BIDS';
% dicomdir = 'E:\mytrainingdata\data exchange server\RABIDS-example\dicomdir';
% subject = 'RaBIDS01';
% suff = '';
% ses_id = 'ses-01';
% study_identifier = 'PSM_BI-STUDIE';
% data_analysis_path = 'E:\mytrainingdata\your project directory';
% series = 2;
% series_n = 't1_mpr_ns_sag_pat2_iso_asy';
% task = 'mprage';
% addsub = 'yes';
% overwrite = 'yes';
%%

data_dir =  'dataset';

dum = 1;

if strcmp(HowExpectDicoms,'allinone')
        dicomd = dicomdir;
elseif strcmp(HowExpectDicoms,'BIDS')
    if strcmp(ses_id,'none')
        dicomd = [dicomdir,filesep,subject];
    else
        dicomd = [dicomdir,filesep,subject,filesep,ses_id];
    end
else
    out{dum,:} = 'User input to object type ''dicom'' not allowed.\nWarning #7\nProgram stops.\n';
    return
end

%% Dicom import

if contains(addsub,'y')
    prefix = 'sub-';
else
    prefix = '';
end

if strcmp(ses_id,'none')
    write_ses = '_';
    subject_dir = [data_analysis_path,filesep,data_dir,filesep,prefix,subject];
else
    write_ses = ['_',ses_id,'_'];
    subject_dir = [data_analysis_path,filesep,data_dir,filesep,prefix,subject,filesep,ses_id];
end

if isfile([subject_dir,filesep,'anat',filesep,prefix,subject,write_ses,'T1w.nii']) || isfile([subject_dir,filesep,'anat',filesep,prefix,subject,write_ses,'T1w.nii.gz'])
    fprintf('Found existing anatomy nifti-file.\n')
    if strcmp(overwrite,'yes')
        fprintf('User permission given to overwrite files.\n');
    else
        fprintf('Permission to overwrite files denied.\n');
        return
    end
end

try

    filter_mprage = [subject,suff,'.MR.',study_identifier,'.',sprintf('%04d',series),'.*'];
    dicm2nii([dicomd,filesep,filter_mprage], [subject_dir,filesep,'anat'], '.nii 3D');
    delete([subject_dir,filesep,'anat',filesep,'dcmHeaders.mat'])

    %% Deface anatomical image and rename to BIDS format
    fprintf('Deface...');
    anat_nii = spm_select('FPList',[subject_dir,filesep,'anat'],[series_n,'.nii']);
    matlabbatch{1}.spm.util.deface.images = {anat_nii};
    spm_jobman('run', matlabbatch);
    delete(anat_nii);

    fprintf('Rename to BIDS format\n');

    anonnii_file = spm_select('FPList',[subject_dir,filesep,'anat'],[series_n,'.*']);
    [~,fn,~] = fileparts(anonnii_file(1,:));
    
    oldf = [subject_dir,filesep,'anat',filesep,fn,'.nii'];
    newf = [subject_dir,filesep,'anat',filesep,prefix,subject,write_ses,'T1w.nii'];    
    movefile(oldf,newf);
    gzip(newf);
    delete(newf);
    
    [~,fn,~] = fileparts(anat_nii(1,:));
    try
        movefile([subject_dir,filesep,'anat',filesep,fn,'.json'],[subject_dir,filesep,'anat',filesep,prefix,subject,write_ses,'T1w.json']);
    catch
        out{dum,:} = 'Nifit-supporting json-file not found.\nSwitch on json-output via dicm2nii.\nWarning #10\n\n';
        dum = dum + 1;
    end

    out{dum,:} = 'Dicom import successful!\n\n';
    return
        
catch
    out{dum,:} = 'Scans for this task were not found for subject/session.\nWarning #8\n';
    return
end