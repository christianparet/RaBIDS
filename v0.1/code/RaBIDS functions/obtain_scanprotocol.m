function [out,scanprotocol] = obtain_scanprotocol(HowExpectDicoms,dicomdir,subject,suff,ses_id,study_identifier,data_analysis_path,max_series,write)
% v0.1 release

% Read headers of MRI images and obtain meta-data information on task and series

% Christian Paret, Central Institute of Mental Health Mannheim, 2018-2020

%% Comment out if in use
% clear
% HowExpectDicoms = 'allinone';
% dicomdir = 'dicomdir';
% subject = 'EFP09';
% suff = '';
% ses_id = 'ses-pre';
% study_identifier = '_';
% data_analysis_path = 'S:\AG_Selfregulation\Projects\EFPTest\Data_analysis';
% max_series = 10;
% write = 'yes';
% 
%% Read data and find volumes of your scan series

if strcmp(HowExpectDicoms,'allinone')
        dicomd = [data_analysis_path,filesep,dicomdir];
elseif strcmp(HowExpectDicoms,'BIDS')
    if strcmp(ses_id,'none')
        dicomd = [dicomdir,filesep,subject];
    else
        dicomd = [dicomdir,filesep,subject,filesep,ses_id];
    end
else
    out{1} = 'User input to object type ''dicom'' not allowed. Must be either ''BIDS'' or ''allinone''.\nProgram stops.\n';
    return
end

k=1;
scanprotocol.name = {'no scan found'};

all_files=dir(strcat(dicomd,filesep,subject,suff,'.MR*.IMA'));
disp('Reading data...');

if isempty(all_files)
    out{k,1} = 'No data found for this session\n\n';   
    return
end

for i=1:max_series
    scan_serie=num2str(i,'%04d');
    ind_files{i}=dir([dicomd,filesep,subject,suff,'.MR.',study_identifier,'.',scan_serie,'.*.IMA']);
end

%% Collect header information

for i=1:length(ind_files)
    series(i,1) = i;
    try
        hdr = dicominfo(fullfile(dicomd,ind_files{i}(1).name));
        if ~strcmp(hdr.PatientID,subject)
            out{k,1} = ['Meta-data ID ''',hdr.PatientID,''' differs from ID assigned in datasheet ''',subject,'''.\nIf this is unexpected you should check whether the volumes are correctly assigned to this subject.\n'];
            k=k+1;
        end
        names(i,1) = {hdr.ProtocolName};
        vols(i,1) = length(ind_files{i});
    catch
        names(i,1) = {'No such series'};
    end
end

%% Save scan protocol to dicomdir

try
    series = series(1:length(vols));
    name = names(1:length(vols));
    scanprotocol = table(series,name,vols);
    if strcmp(write,'yes')
        writetable(scanprotocol,fullfile(dicomd,['scanprotocol_',subject,suff]),'Delimiter','tab')
    end

    out{k,1} = 'Scan protocol saved to dicomdir.\n\n';
    return
    
catch
    out{k,1} = 'Not able to write scan protocol.\n';
    k=k+1;
    
    try
        scanprotocol = readtable(fullfile(dicomd,['scanprotocol_',subject,suff]));
        out{k,1} = 'Found scan protocol in dicomdir, will use that one.\n\n ';
        return
    catch
        out{k,1} = 'No scan protocol found in dicomdir.\nIf this is unexpected, follow steps below:\nIt is recommended to check user input in the datasheet of these object types: data exchange path, dicoms, series info, general suffix and session info.\n';
        return
    end
end