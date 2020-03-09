function out = import2nifti_generic(subject,suff,ses_id,study_identifier,data_analysis_path,max_series,n_images,addsub)
% v0.1 release

% Uses dicm2nii toolbox to import dicom-images to nifti-format
% Folder structure in BIDS format is produced, for more information see: Gorgolewski, K. J. et al., Sci. Data 3:160044, doi: 10.1038/sdata.2016.44 (2016)

% Christian Paret, Central Institute of Mental Health Mannheim, 2018-2020

%% Remove when function in use
% clear
% data_analysis_path = 'S:\AG_Selfregulation\Projects\EFPTest\Data_analysis';
% max_series = 9; % Run N loops; N needs to be equal or larger maximum number of scan series in a subject
% subject=input('Enter subject code (e.g. EFP01)');
% suff = ''; % enter suffix for unambigous identification, e.g. in case two sessions were recorded and were labled with similar name
% ses_id = 'ses-pre'; % enter additional session information for folder structure. if only one session enter 'none'!
% study_identifier = 'EFPNFB';
% n_images = 192; % number of volumes
%%

dicom_dir =  'dicomdir';
data_dir =  'dataset';

%% Read data and find volumes of your scan series

all_files=dir(strcat(data_analysis_path,filesep,dicom_dir,filesep,subject,suff,'.MR*.IMA'));
disp('Reading data...');

if isempty(all_files)
    out = 'No data found for this session\n\n';
    return
end

for i=1:max_series
    scan_serie=num2str(i,'%04d');
    ind_files{i}=dir([data_analysis_path,filesep,dicom_dir,filesep,subject,suff,'.MR.*_',study_identifier,'.',scan_serie,'.*.IMA']);
end
       
% info_file.vp=subject;
% 
% for i=1:max_series
%     if length(ind_files{i})==n_images
%         info_file.mprage=i;
%         break
%     end
% end

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

try

    filter_mprage = [subject,suff,'.MR.*_',study_identifier,'.',sprintf('%04d',info_file.mprage),'.*'];
    dicm2nii([data_analysis_path,filesep,dicom_dir,filesep,filter_mprage], [subject_dir,filesep,'anat'], '.nii 3D');
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
        
catch
    out = 'No such task or series, check label and min/max number of volumes.\n\n';
end