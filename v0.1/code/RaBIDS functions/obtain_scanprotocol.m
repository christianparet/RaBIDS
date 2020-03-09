function [out,scanprotocol] = obtain_scanprotocol(HowExpectDicoms,dicomdir,subject,suff,ses_id,study_identifier,data_analysis_path,max_series,write)
% v0.1 release

% Read headers of MRI images and obtain meta-data information on task and series

% Christian Paret, Central Institute of Mental Health Mannheim, 2018-2020

%% Comment out if in use
% clear
% HowExpectDicoms = 'BIDS';
% dicomdir = 'S:\AG-Austausch\RaBIDS\example\data exchange server\RABIDS-example\dicomdir';
% subject = 'RaBIDS01';
% suff = '';
% ses_id = 'ses-01';
% study_identifier = 'PSM_BI-STUDIE';
% data_analysis_path = 'S:\AG-Austausch\RaBIDS\example\your project directory';
% max_series = 10;
% write = 'yes';

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
    out{k} = 'User input to object type ''dicom'' not allowed. Must be either ''BIDS'' or ''allinone''.\nProgram stops.\n';
    return
end

k=1;

all_files=dir(strcat(dicomd,filesep,subject,suff,'.MR*.IMA'));
disp('Reading data...');

if isempty(all_files)
    out{k} = 'No data found for this session\n\n';   
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
            out{k} = ['Meta-data ID ''',hdr.PatientID,' differs from ID assigned in datasheet ''',subject,'.\n'];
            k=k+1;
        end
        names(i,1) = {hdr.ProtocolName};
        vols(i,1) = length(ind_files{i});
    catch
        names(i,1) = {'No such series'};
    end
end

%% Save scan protocol to dicomdir

series = series(1:length(vols));
name = names(1:length(vols));
scanprotocol = table(series,name,vols);
if strcmp(write,'yes')
    writetable(scanprotocol,fullfile(dicomd,['scanprotocol_',subject,suff]),'Delimiter','tab')
end

out{k} = 'scan protocol ready\n\n';