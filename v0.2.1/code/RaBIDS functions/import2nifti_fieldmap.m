function out = import2nifti_fieldmap(HowExpectDicoms,dicomdir,subject,suff,ses_id,study_identifier,data_analysis_path,series,series_n,TE1,TE2,IntendedFor,addsub,overwrite)
% v0.2.1 release

% Uses dicm2nii toolbox to import dicom-images to nifti-format
% Folder structure in BIDS format is produced, for more information see: Gorgolewski, K. J. et al., Sci. Data 3:160044, doi: 10.1038/sdata.2016.44 (2016)

% CAVE: This script is developed for fieldmap scans from Siemens scanners
% with output of one phase-difference image and two magnitude images. The
% header of the phase-difference image lacks information for TE2.
% Therefore, TE information needs to be entered manually to the datasheet
% file to be added to the json metadata-file. Without adding TE information to
% the json-file, the BIDS validator will return an error message for the
% phase-difference map.

% Christian Paret, Central Institute of Mental Health Mannheim, 2020

% Function is new with v0.2

% Change log to v0.2.1
% - receive TE information for phasedifference image
% - account for zipped nii files nii.gz

%% Comment out if function in use
% HowExpectDicoms = 'BIDS';
% dicomdir = 'E:\mytrainingdata\data exchange server\RABIDS-example\dicomdir';
% subject = 'RaBIDS01';
% suff = '';
% ses_id = 'ses-01';
% study_identifier = 'PSM_BI-STUDIE';
% data_analysis_path = 'E:\mytrainingdata\your project directory';
% series = 4;
% series_n = 'gre_field_mapping';
% TE1 = 4.92;
% TE2 = 7.38;
% IntendedFor{1,1} = 'scenes';
% IntendedFor{1,2} = 'faces';
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
    
    dicm2nii([dicomd,filesep,filter_fm,'*'], [subject_dir,filesep,'fmap'], 'nii.gz');
    delete([subject_dir,filesep,'fmap',filesep,'dcmHeaders.mat'])
    
    nii_file = [subject_dir,filesep,'fmap',filesep,series_n,'*.gz'];  
    get_file = dir(nii_file);
    
    if length(get_file)>1
        out{dum,:} = 'Expected one imported nifti file but found more than one.\nError #13.\nProgram stops.\n';
        return
    end
    
    [~,fn,fe] = fileparts(get_file.name);
    if strcmp(fn,[series_n,'_phase.nii'])
        fieldmaptype = 'phasediff'; 
        if TE1>0 && TE2>0
            out{dum,:} = 'Phase difference-image found. TE specified.\n';
        else
            out{dum,:} = 'Phase difference-image found. Lacking TE information.\nError #11.\n';
        end
        dum = dum + 1;
    else
        fieldmaptype = 'magnitude';
    end
    
    nii = nii_tool('load', [get_file.folder,filesep,get_file.name]);
    ext = strfind(fn,'.nii');
    fname = [subject_dir,filesep,'fmap',filesep,fn(1:ext-1),'.json'];
    
    switch fieldmaptype
        case 'phasediff'
            % Process image
            nii_tool('save', nii, [subject_dir,filesep,'fmap',filesep,prefix,subject,write_ses,'phasediff',fe]) % save volume 1; need to add run-<index> key/value pair or acq-<label> key/value pair to have more than one fmap per sessions
                
            % Adapt information in json file, rename json file
            try
                jsonf = jsondecode(fileread(fname));
                
                % anonymize
                jsonf.PatientName = subject;
                jsonf = rmfield(jsonf,{'PatientSex','PatientAge','PatientSize','PatientWeight','AcquisitionDateTime'});
                
                % add TE
                jsonf = rmfield(jsonf,'EchoTime');
                jsonf.EchoTime1 = TE1;
                jsonf.EchoTime2 = TE2;
                
                % add IntendedFor
                for i = 1:length(IntendedFor)
                    if ischar(IntendedFor{1,i})
                        if strcmp(ses_id,'none')
                            taskp = ['func/',prefix,subject,write_ses,'task-',IntendedFor{1,i},'_bold.nii.gz'];
                        else
                            taskp = [ses_id,'/func/',prefix,subject,write_ses,'task-',IntendedFor{1,i},'_bold.nii.gz'];
                        end
                        jsonf.IntendedFor{i,1} = taskp;
                        out{dum,:} = ['Fieldmap intended for task ',IntendedFor{1,i},'.\n'];
                    else
                        if i == 1
                            out{dum,:} = 'No TaskName assigned to fieldmap phasediff.\nError #12\n';
                        end
                    end
                    dum = dum + 1;
                end
                
                jsonf_newname = [prefix,subject,write_ses,'phasediff.json'];
                saveJSONfile(jsonf,fullfile(subject_dir,filesep,'fmap',filesep,jsonf_newname)) % Lior Kirsch (2020). Structure to JSON (https://www.mathworks.com/matlabcentral/fileexchange/50965-structure-to-json), MATLAB Central File Exchange. Retrieved February 27, 2020.
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
    end   
    
    if isfile(fname) % delete json file; according to BIDS v1.4.1 standard no json file for magnitude image(s) needed
        delete(fname)
    end    
    delete([subject_dir,filesep,'fmap',filesep,fn,'.gz']); % delete non-BIDS conform nii-image
        
    out{dum,:} = 'Dicom import successful!\n\n';
    
catch
    out{dum,:} = 'Scans for this task were not found for this subject/session.\nError #8\n';
    return
end