function out = import2nifti_fieldmap(HowExpectDicoms,dicomdir,subject,suff,ses_id,study_identifier,data_analysis_path,series_counter,series_name,IntendedFor,FieldIdentifier,addsub,addses,overwrite,dcm2niix_exe_path)
% Most recently updated for  v0.4 release

% Uses dicm2nii toolbox to import dicom-images to nifti-format
% Folder structure in BIDS format is produced, for more information see: Gorgolewski, K. J. et al., Sci. Data 3:160044, doi: 10.1038/sdata.2016.44 (2016)

% CAVE: This script is developed for fieldmap scans from Siemens scanners
% with output of one phase-difference image and two magnitude images. The
% header of the phase-difference image lacks information for TE2.
% Therefore, TE information needs to be entered manually to the datasheet
% file. RaBIDS will add this information to the json metadata-file. Without adding TE information to
% the json-file, the BIDS validator will return an error message for the
% phase-difference map.

% Christian Paret, Central Institute of Mental Health Mannheim, 2020

% Function is new with v0.2

% Change log
% from v0.2.2 to v0.4
% - added 'HowExpectDicoms' = 'Siemens' option for new Siemens dicom format
% - option to add session prefix 'ses-'
% - json anonymisation implemented
% - use jsonencode to write json file
% -  a bunch of different export format types of dcms are incorporated (i.e. Siemens_MA, Siemens_TB, Siemens_FR), corresponding to individual idiosyncracies of dcm export parameters at different MR centers.

% v0.2.2
% - remove field TotalReadoutTime from fmap-accompanying json-sidecar file

% v0.2.1
% - receive TE information for phasedifference image
% - account for zipped nii files nii.gz

%% Comment out if function in use
% clear
% HowExpectDicoms = 'Siemens_MA';
% dicomdir = 'W:\group_psm\BrainSTEADy_data_analysis\Neuroimaging\dicoms';
% subject = 'MA101';
% suff = '';
% ses_id = 'V2';
% study_identifier = '';
% data_analysis_path = 'W:\group_psm\BrainSTEADy_data_analysis\Neuroimaging';
% series_counter = 20;
% series_name = 'gre_field_mapping';
% IntendedFor{1,1} = 'task-neurofeedback_run-1';
% IntendedFor{1,2} = 'task-neurofeedback_run-2';
% FieldIdentifier = 'gre_fmap0';
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
    dum = dum + 1;
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
    write_ses2 = '/';
else
    if contains(addses,'y')
        prefix_ses = 'ses-';
    else
        prefix_ses = '';
    end
    write_ses = ['_',prefix_ses,ses_id,'_'];
    write_ses2 = ['/',prefix_ses,ses_id,'/'];
    subject_dir = [data_analysis_path,filesep,data_dir,filesep,prefix_sub,subject,filesep,prefix_ses,ses_id];
end

fmapp = fullfile(subject_dir,'fmap');

if ~isfolder(fmapp)
    mkdir(fmapp)
elseif isfile([fmapp,filesep,prefix_sub,subject,write_ses,'phasediff.nii']) || isfile([fmapp,filesep,prefix_sub,subject,write_ses,'phasediff.nii.gz'])
    fprintf('Found existing fieldmap nifti-file.\n')
    if strcmp(overwrite,'yes')
        % Need to delete existing files
        fprintf('User permission given to overwrite files.\n');
        dum1 = dir(fullfile(fmapp,'sub-*'));
        nr_niftis = length(dum1);
        for i = 1:nr_niftis
            delete(fullfile(fmapp,dum1(i).name))
        end
    else
         out{dum,:} = 'Permission to overwrite files denied.\n\n';
        return
    end
end

try
    if strcmp(HowExpectDicoms,'Siemens_MA')
        filter_fm = [series_name,'_',num2str(series_counter),'_MR'];
    elseif strcmp(HowExpectDicoms,'Siemens_TB')
        filter_fm = [sprintf('%02d',series_counter),'_',series_name];
    elseif strcmp(HowExpectDicoms,'Siemens_FR')
        filter_fm = [sprintf('%d',series_counter),'_',lower(series_name)];
    else
        filter_fm = [subject,suff,'.MR.',study_identifier,'.',sprintf('%04d',series_counter),'.*'];
    end

    dcm_input = fullfile(dicomd,filter_fm);

    dum1 = dir(fullfile(dcm_input,'*'));
    nr_dcms = length(dum1(not([dum1.isdir])));

    if nr_dcms == 1 % Workaround: phasediff is single image, whereas two magnitude images are expected
        nifti_name = [prefix_sub,subject,write_ses,'phasediff'];
    elseif nr_dcms == 2
        nifti_name = [prefix_sub,subject,write_ses,'magnitude'];
    else
        out{dum,:} = 'There are more fieldmap dicom images as expected.\n\n';
        dum = dum + 1;
        return
    end

    system([dcm2niix_exe_path ... % path to executable
        ' -a "y"' ... % adjacent DICOMs (images from same series always in same folder) for faster conversion (n/y, default n)
        ' -f ' , nifti_name , ... % filename (%a=antenna (coil) name, %b=basename, %c=comments, %d=description, %e=echo number, %f=folder name, %g=accession number, %i=ID of patient, %j=seriesInstanceUID, %k=studyInstanceUID, %m=manufacturer, %n=name of patient, %o=mediaObjectInstanceUID, %p=protocol, %r=instance number, %s=series number, %t=time, %u=acquisition number, %v=vendor, %x=study ID; %z=sequence name; default '%f_%p_%t_%s')
        ' -o ' , fmapp , ... % output directory
        ' -z "y" ' ... % -z : gz compress images (y/i/n/3, default n) [y=pigz, i=internal:zlib, n=no, 3=no,3D]
        dcm_input]); % location of dicoms

    niftis = dir(fullfile(fmapp,[nifti_name,'*']));

    % Naming convention of dcm2niix fieldmap output is not BIDS valid. We
    % are using a workaround to rename files, based on assumption, that 1
    % phasedifference image and 2 magnitude images are available
    if nr_dcms == 1
        pos = strfind(niftis(1).name,'phasediff');
        for i = 1:length(niftis)
            new_nifti_name = [niftis(i).name(1:pos+8),niftis(i).name(pos+15:end)];
            movefile(fullfile(niftis(i).folder,niftis(i).name), fullfile(niftis(i).folder,new_nifti_name));
        end
    elseif nr_dcms == 2
        pos = strfind(niftis(1).name,'magnitude');
        for i = 1:length(niftis)
            new_nifti_name = [niftis(i).name(1:pos+8),niftis(i).name(pos+11:end)];
            movefile(fullfile(niftis(i).folder,niftis(i).name), fullfile(niftis(i).folder,new_nifti_name));
        end
    end

    out{dum,:} = 'Dicom import successful!\n\n';

    if nr_dcms == 1

        json_info = jsondecode(fileread(fullfile(fmapp,[nifti_name,'.json'])));
    
        if any(FieldIdentifier)
            json_info.B0FieldIdentifier = [FieldIdentifier,write_ses];
            out{dum,:} = ['B0FieldIdentifier ',[FieldIdentifier,write_ses],' added to json-file.\n'];
            dum = dum + 1;
        else
            out{dum,:} = 'B0FieldIdentifier not specified.\n';
            dum = dum + 1;
        end
    
        for i = 1:length(IntendedFor)
            if ischar(IntendedFor{1,i})
                taskp = ['bids::',prefix_sub,subject,write_ses2,'func/',prefix_sub,subject,write_ses,IntendedFor{1,i},'_bold.nii.gz'];
                json_info.IntendedFor{i,1} = taskp;
                out{dum,:} = ['Fieldmap intended for ',IntendedFor{1,i},'.\n'];
            else
                if i == 1
                    out{dum,:} = 'No TaskName assigned to fieldmap phasediff.\nWarning #12\n';
                end
            end
            dum = dum + 1;
        end
    
        json_new = jsonencode(json_info,PrettyPrint=true);
        fid = fopen(fullfile(fmapp,[nifti_name,'.json']), 'w');
        fprintf(fid, '%s', json_new);
        fclose(fid);
    
    end

catch
    out{dum,:} = 'Scans for this task were not found for this subject/session.\nWarning #8\n';
    dum = dum + 1;
    return
end