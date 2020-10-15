function out = import2nifti_fieldmap(HowExpectDicoms,dicomdir,subject,suff,ses_id,study_identifier,data_analysis_path,series,series_n,addsub,overwrite)
% v0.2 release

% Uses dicm2nii toolbox to import dicom-images to nifti-format
% Folder structure in BIDS format is produced, for more information see: Gorgolewski, K. J. et al., Sci. Data 3:160044, doi: 10.1038/sdata.2016.44 (2016)

% Current version only supports one fieldmap per session
%% Supports fieldmaps with phasedifference image only. DOUBLE-CHECK WHETHER IT IS REALLY A PHASEDIFFERENCE I AM USING!!!

% Christian Paret, Central Institute of Mental Health Mannheim, 2020

% New with v0.2

%% Comment out if function in use
% HowExpectDicoms = 'BIDS';
% dicomdir = 'E:\mytrainingdata\data exchange server\RABIDS-example\dicomdir';
% subject = 'RaBIDS01';
% suff = '';
% ses_id = 'ses-01';
% study_identifier = 'PSM_BI-STUDIE';
% data_analysis_path = 'E:\mytrainingdata\your project directory';
% series = 3;
% series_n = 'gre_field_mapping';
% addsub = 'yes';
% overwrite = 'yes';
%%

data_dir =  'dataset';

dum = 1;

if strcmp(HowExpectDicoms,'allinone')
        dicomd = [data_analysis_path,filesep,dicomdir];
elseif strcmp(HowExpectDicoms,'BIDS')
    if strcmp(ses_id,'none')
        dicomd = [dicomdir,filesep,subject];
    else
        dicomd = [dicomdir,filesep,subject,filesep,ses_id];
    end
else
    out{dum,:} = 'User input to object type ''dicom'' not allowed. Must be either ''BIDS'' or ''allinone''.\nProgram stops.\n';
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

if isfolder([subject_dir,filesep,'fmap'])
    fprintf('Fieldmap directory exists.\n')
    if strcmp(overwrite,'yes')
        fprintf('User permission given to write into fmap directory.\n');
    else
        fprintf('Permission to overwrite files denied.\n');
        return
    end
end

try
    
    filter_fm = [subject,suff,'.MR.',study_identifier,'.',sprintf('%04d',series),'.'];
    
    dicm2nii([dicomd,filesep,filter_fm,'*'], [subject_dir,filesep,'fmap'], '.nii');
    delete([subject_dir,filesep,'fmap',filesep,'dcmHeaders.mat'])
    
    nii_file = [subject_dir,filesep,'fmap',filesep,series_n,'*.nii'];  
    get_file = dir(nii_file);
    [~,fn,fe] = fileparts(get_file.name);
    if strcmp(fn,[series_n,'_phase'])
        fieldmaptype = 'phasediff'; % Is this really a phase diff image???
        out{dum,:} = 'Phase-image found.\n';
        dum = dum + 1;
    else
        fieldmaptype = 'magnitude';
    end
    
    nii = nii_tool('load', [get_file.folder,filesep,get_file.name]);
    
    switch fieldmaptype
        case 'phasediff'
            % Process image
            nii_tool('save', nii, [subject_dir,filesep,'fmap',filesep,prefix,subject,write_ses,'phasediff',fe]) % save volume 1; need to add run-<index> key/value pair or acq-<label> key/value pair to have more than one fmap per sessions

            % Rename json file
            try
                movefile([subject_dir,filesep,'fmap',filesep,fn,'.json'],[subject_dir,filesep,'fmap',filesep,prefix,subject,write_ses,'phasediff.json']);
            catch
                out{dum,:} = 'json file for phasedifference map not found. Consider checking dicm2nii json-output options.\nError #10\n\n';
                dum = dum + 1;
            end

        case 'magnitude'
            % Process image(s)
            nii1 = nii;
            nii1.img = nii.img(:,:,:,1); % take first volume
            nii_tool('save', nii1, [subject_dir,filesep,'fmap',filesep,prefix,subject,write_ses,'magnitude1',fe]) % save volume 1; need to add run-<index> key/value pair or acq-<label> key/value pair to have more than one fmap per sessions

            try 
                nii.img = nii.img(:,:,:,2); % take second volume
                out{dum,:} = 'Two magnitude images found.\n';
                nii_tool('save', nii, [subject_dir,filesep,'fmap',filesep,prefix,subject,write_ses,'magnitude2',fe]) % save volume 2; need to add run-<index> key/value pair or acq-<label> key/value pair to have more than one fmap per sessions  
            catch
                out{dum,:} = 'Found a single magnitude image.\n';
            end
            dum = dum + 1;
            
            if isfile([subject_dir,filesep,'fmap',filesep,fn,'.json']) % delete json file; according to BIDS v1.4.1 standard no json file for magnitude image(s) needed
                delete([subject_dir,filesep,'fmap',filesep,fn,'.json'])
            end
                
    end   
    
    delete([subject_dir,filesep,'fmap',filesep,fn,'.nii']); % delete non-BIDS conform nii-image
        
    out{dum,:} = 'Dicom import successful!\n\n';
    
catch
    out{dum,:} = 'Scans for this task were not found for this subject/session.\nError #8\n';
    return
end