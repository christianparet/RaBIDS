function [out,scanprotocol] = obtain_scanprotocol(HowExpectDicoms,dicomdir,subject,suff,ses_id,study_identifier,data_analysis_path,max_series,addsub,addses,write)

% Most recently updated for v0.4 release

% Read headers of MRI images and obtain meta-data information on task and series

% Christian Paret, Central Institute of Mental Health Mannheim, 2018-2025

% Change log from v0.1:
% - writes to dataset directory
% - sanity check of MR series names
% - create subject-directory if not existant

% Change log from v0.2 to v0.4
% - added 'HowExpectDicoms' = 'Siemens' option for new Siemens dicom format
% - option to add session prefix 'ses-'
% - a bunch of different export format types of dcms are incorporated (i.e. Siemens_MA, Siemens_TB, Siemens_FR), corresponding to individual idiosyncracies of dcm export parameters at different MR centers.

%% Comment out if in use
% clear
% HowExpectDicoms = 'Siemens_FR';
% dicomdir = 'W:\group_psm\BrainSTEADy-Datenanalyse\MRT\dicoms';
% subject = 'FR103';
% suff = '';
% ses_id = 'V2';
% study_identifier = '';
% data_analysis_path = 'W:\group_psm\BrainSTEADy-Datenanalyse\MRT';
% max_series = 30;
% addsub = 'yes';
% addses = 'yes';
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

%% Read data and find volumes of your scan series

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

scanprotocol.name = {'no scan found'};

if strcmp(HowExpectDicoms,'Siemens_MA')

    for i=1:max_series
        scan_serie=num2str(i);
        ind_files{i}=dir(fullfile(dicomd,['*_',scan_serie,'_MR']));
    end

elseif strcmp(HowExpectDicoms,'Siemens_TB')

    for i=1:max_series
        ind_files{i}=dir(fullfile(dicomd,[sprintf('%02d',i),'_*']));
    end

elseif strcmp(HowExpectDicoms,'Siemens_FR')

    for i=1:max_series
        ind_files{i}=dir(fullfile(dicomd,[sprintf('%d',i),'_*']));
    end

else

    all_files=dir(strcat(dicomd,filesep,subject,suff,'.MR*.IMA'));
    disp('Reading data...');

    if isempty(all_files)
        out{dum,1} = 'No data found for this session\n.Warning #12\n\n';
        return
    end

    for i=1:max_series
        scan_serie=num2str(i,'%04d');
        ind_files{i}=dir([dicomd,filesep,subject,suff,'.MR.',study_identifier,'.',scan_serie,'.*.IMA']);
    end

end

%% Collect header information
scaninfofile = fullfile(subject_dir,[prefix_sub,subject,write_ses,'ScanInfo.mat']);
if exist(scaninfofile) == 2
    fprintf('Found existing ScanInfo file including header information. Will use this file.\n')
    load(scaninfofile)
else

    fprintf('Reading image headers ...\n'); % To do: implement routine to save and load hdr file
    
    for i=1:length(ind_files)
        series(i,1) = i;

        try
            if strcmp(HowExpectDicoms,'Siemens_MA')
                hdr = dicominfo(fullfile(dicomd,ind_files{i}(1).name,'1.dcm'));
            elseif strcmp(HowExpectDicoms,'Siemens_TB')
                if isfile(fullfile(dicomd,ind_files{i}(1).name,'0001_0001'))
                    hdr = dicominfo(fullfile(dicomd,ind_files{i}(1).name,'0001_0001'));
                elseif isfile(fullfile(dicomd,ind_files{i}(1).name,'0_0001'))
                    hdr = dicominfo(fullfile(dicomd,ind_files{i}(1).name,'0_0001'));
                end
            elseif strcmp(HowExpectDicoms,'Siemens_FR')
                dcm = dir(fullfile(dicomd,ind_files{i}(1).name,'MRe.*.dcm'));
                hdr = dicominfo(fullfile(dicomd,ind_files{i}(1).name,dcm(1).name));
            else
                hdr = dicominfo(fullfile(dicomd,ind_files{i}(1).name));
            end
            if ~strcmp(hdr.PatientID,subject)
                out{dum,1} = ['MRI-series #',num2str(i),': meta-data ID ''',hdr.PatientID,''' and ID assigned in datasheet ''',subject,''' are different.\nWarning #1.\n\n'];
                dum=dum+1;
            end
            names(i,1) = {hdr.ProtocolName};
    
            if strcmp(HowExpectDicoms,'Siemens_MA') || strcmp(HowExpectDicoms,'Siemens_FR')
                dcmfiles = dir(fullfile(dicomd,ind_files{i}(1).name,'*.dcm'));
                vols(i,1) = length(dcmfiles);
            elseif strcmp(HowExpectDicoms,'Siemens_TB')
                dcmfiles = dir(fullfile(dicomd,ind_files{i}(1).name,'0*'));
                vols(i,1) = length(dcmfiles);
            else
                vols(i,1) = length(ind_files{i});
            end
        catch
            names(i,1) = {'No such series'};
        end
    end

    %% Safe header and volume information
    if strcmp(write,'yes')
        if ~isfolder(subject_dir)
            mkdir(subject_dir)
        end
        save(scaninfofile,"hdr","names","vols","series");
        fprintf('Scan information file ScanInfo was saved to subject directory.\n');
    end
end

%% Sanity check
fprintf('Checking for duplicate names of sequences ...\n');

[D,~,IC] = unique(names(:));
Y = hist(IC,unique(IC));
Z = struct('name',D,'freq',num2cell(Y(:)));

[unique_seq,~] = size(Z);

for i=1:unique_seq
    if ~strcmp(Z(i).name,'No such series')
        if Z(i).freq>1
            out{dum,1} = ['Found ',num2str(Z(i).freq),' MR series with name ',Z(i).name,'.\nFor BOLD series with identical names: ranges defined by MinImages-MaxImages must not overlap.\nWarning#11\n\n'];
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
        writetable(scanprotocol,fullfile(subject_dir,[prefix_sub,subject,write_ses,'scanprotocol.txt']),'Delimiter','tab')
    end

    out{dum,1} = 'Scan protocol saved.\n\n';
    return
    
catch
    out{dum,1} = 'Not able to save scan protocol.\n';
    dum=dum+1;

    try
        scanprotocol = readtable(fullfile(subject_dir,[prefix_sub,subject,write_ses,'scanprotocol.txt']));
        out{dum,1} = 'Found scanprotocol, will use that one.\n\n ';
        return
    catch
        out{dum,1} = 'No scan protocol found.\nWarning #2.\n\n';
        return
    end
end