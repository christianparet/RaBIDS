function out = import2nifti_anatomical(HowExpectDicoms,dicomdir,subject,suff,ses_id,study_identifier,data_analysis_path,series,addsub,overwrite)
% v0.1 release

% Uses dicm2nii toolbox to import dicom-images to nifti-format
% Folder structure in BIDS format is produced, for more information see: Gorgolewski, K. J. et al., Sci. Data 3:160044, doi: 10.1038/sdata.2016.44 (2016)

% Christian Paret, Central Institute of Mental Health Mannheim, 2018-2020


%%

data_dir =  'dataset';

if strcmp(HowExpectDicoms,'allinone')
        dicomd = [data_analysis_path,filesep,dicomdir];
elseif strcmp(HowExpectDicoms,'BIDS')
    if strcmp(ses_id,'none')
        dicomd = [dicomdir,filesep,subject];
    else
        dicomd = [dicomdir,filesep,subject,filesep,ses_id];
    end
else
    out = 'User input to object type ''dicom'' not allowed. Must be either ''BIDS'' or ''allinone''.\nProgram stops.\n';
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

if isfolder([subject_dir,filesep,'anat'])
    fprintf('Directory with functional images exists.\n')
    if strcmp(overwrite,'yes')
        fprintf('User permission given to overwrite.\n');
    else
        fprintf('Permission to overwrite declined.\n');
        return
    end
end

try

    filter_mprage = [subject,suff,'.MR.',study_identifier,'.',sprintf('%04d',series),'.*'];
    dicm2nii([dicomd,filesep,filter_mprage], [subject_dir,filesep,'anat'], '.nii 3D');
    delete([subject_dir,filesep,'anat',filesep,'dcmHeaders.mat'])

    %% Deface anatomical image and rename to BIDS format
    fprintf('Deface\n');
    anat_nii = spm_select('FPList',[subject_dir,filesep,'anat'],'^t1_mpr.*.nii');
    matlabbatch{1}.spm.util.deface.images = {anat_nii};
    spm_jobman('run', matlabbatch);
    delete(anat_nii);

    fprintf('Rename to BIDS format\n');
    nii_file = spm_select('FPList',[subject_dir,filesep,'anat'],'t1_mpr.*');

    [nr_files,~] = size(nii_file); % there will be two files if option to produce .json file is activated via dicm2nii GUI

    for i = 1:nr_files
        [~,fn,fe] = fileparts(nii_file(i,:));
        movefile([subject_dir,filesep,'anat',filesep,fn,fe],[subject_dir,filesep,'anat',filesep,prefix,subject,write_ses,'T1w',fe]);
    end

    out='Dicom import successful!\n\n';
    return
        
catch
    out = 'No such task or series.\nIf this is unexpected, follow steps below:\nIf program was not able to save scan protocol with current settings, it is recommended to check user input in the datasheet of these object types: data exchange path, dicoms, series info, general suffix and session info.\nIf scan protocol of this subject exists in dicomdir, go to datasheet and check the sobject type MRI series, minimum and maximum images, and object type series info.\nCompare user input with scan protocol of this subject, which you find in the dicomdir.\nIs MRI series ID (= name in scan protocol) and number of images (vols in scan protocol) appropriately defined?\n\n';
    return
end