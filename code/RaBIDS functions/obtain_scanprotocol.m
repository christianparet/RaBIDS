function [out,scanprotocol] = obtain_scanprotocol(HowExpectDicoms,dicomdir,subject,suff,ses_id,study_identifier,data_analysis_path,max_series,addsub,write)
% v0.2 release

% Read headers of MRI images and obtain meta-data information on task and series

% Christian Paret, Central Institute of Mental Health Mannheim, 2018-2020

% Change log from v0.1:
% - writes to dataset directory
% - sanity check of MR series names
% - create subject-directory if not existant

%% Comment out if in use
% clear
% HowExpectDicoms = 'BIDS';
% dicomdir = 'E:\mytrainingdata\data exchange server\RABIDS-example\dicomdir';
% subject = 'RABIDS01';
% suff = '';
% ses_id = 'ses-01';
% study_identifier = 'PSM_BI-STUDIE';
% data_analysis_path = 'E:\mytrainingdata\your project directory';
% max_series = 10;
% addsub = 'yes';
% write = 'yes';
%%

data_dir =  'dataset';

dum=1;

if strcmp(HowExpectDicoms,'allinone')
        dicomd = dicomdir;
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

%% Read data and find volumes of your scan series

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

scanprotocol.name = {'no scan found'};

all_files=dir(strcat(dicomd,filesep,subject,suff,'.MR*.IMA'));
disp('Reading data...');

if isempty(all_files)
    out{dum,1} = 'No data found for this session\n.Error #12\n\n';   
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
            out{dum,1} = ['MRI-series #',num2str(i),': meta-data ID ''',hdr.PatientID,''' and ID assigned in datasheet ''',subject,''' are different.\nError #1.\n\n'];
            dum=dum+1;
        end
        names(i,1) = {hdr.ProtocolName};
        vols(i,1) = length(ind_files{i});
    catch
        names(i,1) = {'No such series'};
    end
end

%% Sanity check

[D,~,IC] = unique(names(:));
Y = hist(IC,unique(IC));
Z = struct('name',D,'freq',num2cell(Y(:)));

[unique_seq,~] = size(Z);

for i=1:unique_seq
    if ~strcmp(Z(i).name,'No such series')
        if Z(i).freq>1
            out{dum,1} = ['Found ',num2str(Z(i).freq),' MR series with name ',Z(i).name,'.\nFor BOLD series with identical names: ranges defined by MinImages-MaxImages must not overlap.\nError#11\n\n'];
            dum=dum+1;
        end
    end
end 

%% Save scan protocol to dicomdir

try
    series = series(1:length(vols));
    name = names(1:length(vols));
    scanprotocol = table(series,name,vols);
    if strcmp(write,'yes')
        if ~isfolder(subject_dir)
            mkdir(subject_dir)
        end
        writetable(scanprotocol,fullfile(subject_dir,[prefix,subject,'_',ses_id,'_scanprotocol.txt']),'Delimiter','tab')
    end

    out{dum,1} = 'Scan protocol saved.\n\n';
    return
    
catch
    out{dum,1} = 'Not able to write scan protocol.\n';
    dum=dum+1;
    
    try
        scanprotocol = readtable(fullfile(subject_dir,[prefix,subject,'_',ses_id,'_scanprotocol.txt']));
        out{dum,1} = 'Found scanprotocol, will use that one.\n\n ';
        return
    catch
        out{dum,1} = 'No scan protocol found.\nError #2.\n\n';
        return
    end
end