function out = import2nifti_diffusion(HowExpectDicoms,dicomdir,subject,suff,ses_id,study_identifier,data_analysis_path,series_counter,series_name,acq,direction,taskrun,addsub,addses,overwrite,dcm2niix_exe_path)

%% Comment out if function in use
% clear
% HowExpectDicoms = 'Siemens_MA';
% dicomdir = 'W:\group_psm\BrainSTEADy_data_analysis\Neuroimaging\dicoms';
% subject = 'MA101';
% suff = '';
% ses_id = 'V2';
% study_identifier = '';
% data_analysis_path = 'W:\group_psm\BrainSTEADy_data_analysis\Neuroimaging';
% series_counter = 8;
% series_name = 'dwi_ses-base_acq-cmrr_task-rest_run-02_dir-PA';
% acq = 'cmrr';
% direction = 'PA';
% taskrun = 2;
% addsub = 'yes';
% addses = 'yes';
% overwrite = 'yes';
% dcm2niix_exe_path = 'E:\dcm2niix\dcm2niix.exe';
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
elseif strcmp(HowExpectDicoms,'Siemens_MA') || strcmp(HowExpectDicoms,'Siemens_TB')
    dicomd_temp = fullfile(dicomdir,[subject,'-',ses_id]);
    dumdir1 = dir(dicomd_temp);
    dumdir2 = dir(fullfile(dicomd_temp,dumdir1(3).name));
    dicomd = fullfile(dicomd_temp,dumdir1(3).name,dumdir2(3).name);
elseif strcmp(HowExpectDicoms,'Siemens_FR')
    dicomd_temp = fullfile(dicomdir,[subject,'-',ses_id]);
    dumdir1 = dir(dicomd_temp);
    dicomd = fullfile(dicomd_temp,dumdir1(3).name);
else
    out{dum,:} = 'User input to object type ''dicom'' not allowed.\nWarning #7\nProgram stops.\n';
    return
end

%% Dicom import

if contains(addsub,'y')
    prefix_sub = 'sub-';
else
    prefix_sub = '';
end

if strcmp(ses_id,'none')
    subject_dir = [data_analysis_path,filesep,data_dir,filesep,prefix_sub,subject];
    write_ses = '_';
else
    if contains(addses,'y')
        prefix_ses = 'ses-';
    else
        prefix_ses = '';
    end
    write_ses = ['_',prefix_ses,ses_id,'_'];
    subject_dir = [data_analysis_path,filesep,data_dir,filesep,prefix_sub,subject,filesep,prefix_ses,ses_id];
end

diffp = fullfile(subject_dir,'dwi');
nifti_name = [prefix_sub,subject,write_ses,'acq-',acq,'_dir-',direction,'_run-',sprintf('%01d',taskrun),'_dwi'];

if ~isfolder(diffp)
    mkdir(diffp)
elseif isfile([diffp,filesep,nifti_name,'.nii']) || isfile([diffp,filesep,nifti_name,'.nii.gz'])
    fprintf('Found existing nifti-file for this task.\n');
    if strcmp(overwrite,'yes')
        % Need to delete existing files
        fprintf('User permission given to overwrite files.\n');
        dum1 = dir(fullfile(diffp,'sub-*'));
        nr_niftis = length(dum1);
        for i = 1:nr_niftis
            delete(fullfile(diffp,dum1(i).name))
        end
    else
        out{dum,:} = ('Permission to overwrite files denied.\n\n');
        return
    end
end

try
    
    if strcmp(HowExpectDicoms,'Siemens_MA')
        filter_dwi = [series_name,'_',num2str(series_counter),'_MR'];
    elseif strcmp(HowExpectDicoms,'Siemens_TB')
        filter_dwi = [sprintf('%02d',series_counter),'_',series_name];
    elseif strcmp(HowExpectDicoms,'Siemens_FR')
        filter_dwi = [sprintf('%d',series_counter),'_',lower(series_name)];
    else
        filter_dwi = [subject,suff,'.MR.',study_identifier,'.',sprintf('%04d',series_counter),'.*'];
    end
        
    dcm_input = fullfile(dicomd,filter_dwi);

    system([dcm2niix_exe_path ... % path to executable
        ' -a "y"' ... % adjacent DICOMs (images from same series always in same folder) for faster conversion (n/y, default n) 
        ' -f ' , nifti_name , ... % filename (%a=antenna (coil) name, %b=basename, %c=comments, %d=description, %e=echo number, %f=folder name, %g=accession number, %i=ID of patient, %j=seriesInstanceUID, %k=studyInstanceUID, %m=manufacturer, %n=name of patient, %o=mediaObjectInstanceUID, %p=protocol, %r=instance number, %s=series number, %t=time, %u=acquisition number, %v=vendor, %x=study ID; %z=sequence name; default '%f_%p_%t_%s') 
        ' -o ' , diffp , ... % output directory
        ' -z "y" ' ... % -z : gz compress images (y/i/n/3, default n) [y=pigz, i=internal:zlib, n=no, 3=no,3D] 
        dcm_input]); % location of dicoms
                  
    out{dum,:} = 'Dicom import successful!\n\n';
    
catch
    out{dum,:} = 'Scans for dwi were not found for this subject/session.\nWarning #8\n';
    return
end