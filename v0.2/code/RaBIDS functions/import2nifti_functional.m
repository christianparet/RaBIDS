function out = import2nifti_functional(HowExpectDicoms,dicomdir,subject,suff,ses_id,study_identifier,first_image,data_analysis_path,series,series_n,task,addsub,overwrite)
% v0.2 release

% Uses dicm2nii toolbox to import dicom-images to nifti-format
% Folder structure in BIDS format is produced, for more information see: Gorgolewski, K. J. et al., Sci. Data 3:160044, doi: 10.1038/sdata.2016.44 (2016)

% Christian Paret, Central Institute of Mental Health Mannheim, 2019-2020

% Change log from v0.1:
% - implemented json-output check
% - implemented Error-reporting

%% Comment out if function in use
% HowExpectDicoms = 'BIDS';
% dicomdir = 'E:\mytrainingdata\data exchange server\RABIDS-example\dicomdir';
% subject = 'RaBIDS01';
% suff = '';
% ses_id = 'ses-01';
% study_identifier = 'PSM_BI-STUDIE';
% data_analysis_path = 'E:\mytrainingdata\your project directory';
% series = 6;
% series_n = 'ep2d_TR2000_64_Scenes';
% task = 'scenes';
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
    out{dum,:} = 'User input to object type ''dicom'' not allowed.\nError #7\nProgram stops.\n';
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

if isfile([subject_dir,filesep,'func',filesep,prefix,subject,write_ses,'task-',task,'_bold.nii'])
    fprintf('Found existing nifti-file for this task.\n')
    if strcmp(overwrite,'yes')
        fprintf('User permission given to overwrite files.\n');
    else
        fprintf('Permission to overwrite files denied.\n');
        return
    end
end

try
    
    filter_epi = [subject,suff,'.MR.',study_identifier,'.',sprintf('%04d',series),'.'];
    
    %  delete initial functional images to eliminate T1-effects
    imgcount=0;
    for i=1:first_image-1
        searchfile = dir([dicomd,filesep,filter_epi,sprintf('%04d',i),'.*']);
        fileex = exist([searchfile.folder,filesep,searchfile.name],'file');
        if fileex==2
            imgcount=imgcount+1;
            datastruct{imgcount}=searchfile;
        end
    end
    
    if imgcount
        fprintf('Initial epi-images are being deleted:\n');
        for i=1:length(datastruct)
            fprintf([datastruct{i}.name,'\n']);
        end
        for i=1:length(datastruct)
            try
                delete([datastruct{i}.folder,filesep,datastruct{i}.name]);
            catch
                out{dum,:} = ['Cannot delete ',datastruct{i}.name,', check permission.\nProgram stops\n'];
                return
            end
        end
    else
        fprintf('No initial functional images found, continue.\n');
    end
    
    dicm2nii([dicomd,filesep,filter_epi,'*'], [subject_dir,filesep,'func'], '.nii');
    delete([subject_dir,filesep,'func',filesep,'dcmHeaders.mat'])
    
    %% Rename to BIDS format
    
    fprintf('Rename to BIDS format\n');
       
    nii_file = spm_select('FPList',[subject_dir,filesep,'func'],[series_n,'.*']);
    
    [nr_files,~] = size(nii_file); % there will be two files if option to produce .json file is activated via dicm2nii GUI
    if ~nr_files>1
        out{dum,:} = 'Nifit-supporting json-file not found.\nSwitch on json-output via dicm2nii before you continue.\nError #10\n';
        dum = dum + 1;
    end
    
    for i = 1:nr_files
        [~,fn,fe] = fileparts(nii_file(i,:));
        movefile([subject_dir,filesep,'func',filesep,fn,fe],[subject_dir,filesep,'func',filesep,prefix,subject,write_ses,'task-',task,'_bold',fe]);
    end
    
    out{dum,:} = 'Dicom import successful!\n\n';
    
    % Write task information to json file. This code could be extended in future with more task information provided in the datasheet table
    % Lior Kirsch (2020). Structure to JSON (https://www.mathworks.com/matlabcentral/fileexchange/50965-structure-to-json), MATLAB Central File Exchange. Retrieved February 27, 2020.
    if ~strcmp(task,'rest') 
        jsonfile_name = ['task-',task,'_bold.json'];
        if ~isfile(fullfile(data_analysis_path,data_dir,jsonfile_name))
            jsonf.TaskName = ['task-',task];
            saveJSONfile(jsonf,fullfile(data_analysis_path,data_dir,jsonfile_name))
        end
    end
    
catch
    out{dum,:} = 'Scans for this task were not found for this subject/session.\nError #8\n';
    return
end