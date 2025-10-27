%% RaBIDS - Rapid preparation of BIDS
% Most recently updated for v0.4 release

% Change log 
% Change log to v0.4
% - add session prefix
% - deface option removed
% - opts input to readtable function added to preserve character variable type when reading data from datasheet (Miroslava Jindrova, 2022/03/04)
% - a bunch of different export format types of dcms are incorporated (i.e. Siemens_MA, Siemens_TB, Siemens_FR), corresponding to individual idiosyncracies of dcm export parameters at different MR centers.

% v0.2.2 
% - added fieldmap import

% Change log to v0.2.1
% - for fieldmap phasedifference maps: read TEs from datasheet and push them to fieldmap-nifti-import function

% v0.1: 
% - see obtain_scanprotocol; scan protocol saved to subject directory
% - save command window output to RaBIDS-logs directory
% - new entries added to bidsignore-file to allow above output in BIDS validation
% - Datasheet variable name naming convention changed
% - Implemented fieldmap import
% - Warning reporting with reference to RaBIDS Warning Reference-sheet
% - Removed several bugs


%% Import imaging data to BIDS format
% Program automatically imports your MRI data to nifti and rearranges it according to the BIDS data format

% Christian Paret, ZI Mannheim, 2019-2020; christian.paret@zi-mannheim.de

clear
clc
fprintf('RaBIDS ---------- Rapid preparation of BIDS ------- is running!\nVersion 0.4\n');

%% Read data from datasheet
opts = detectImportOptions('datasheet_MA.xlsx','NumHeaderLines',0);
opts.PreserveVariableNames = 1;
opts = setvartype(opts,3,'char');
data = readtable('datasheet_MA.xlsx',opts,'ReadRowNames',true);
data_dir =  'dataset';

max_series = 70; % N needs to be equal or larger maximum number of scan series in a subject; this could be improved in future versions

descriptcol = find(strcmp(data.Properties.VariableNames,'Description'));
userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));
runNumbercol = find(strcmp(data.Properties.VariableNames,'RunNumber'));
minImagescol = find(strcmp(data.Properties.VariableNames,'MinImages'));
maxImagescol = find(strcmp(data.Properties.VariableNames,'MaxImages'));

AddPathline = find(contains(data.Properties.RowNames,'add path'));
for i = 1:length(AddPathline)
    addpath(data{AddPathline(i),userInputcol}{:});
end

Dcm2niixExePathline = find(contains(data.Properties.RowNames,'executable'));
dcm2niix_exe_path = data{Dcm2niixExePathline,userInputcol}{:};

DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};

Dicomsline = find(strcmp(data.Properties.RowNames,'dicoms'));
HowExpectDicoms = data{Dicomsline,userInputcol}{:};
ExchangeServerPathline = find(strcmp(data.Properties.RowNames,'data exchange path'));
dicomdir = data{ExchangeServerPathline,userInputcol}{:};

dum = find(contains(data.Properties.RowNames,'subject info'));
for i = 1:length(dum)
    subj_list(i).name = data{dum(i),userInputcol};
end

general_suffixline = find(strcmp(data.Properties.RowNames,'general suffix'));
general_suffix = data{general_suffixline,userInputcol}{:}; % study code at scanner site, used to filter dicoms

addsubline = find(strcmp(data.Properties.RowNames,'add sub prefix')); % if subject code does not begin with 'sub-' this needs to be prepended for BIDS compatibility
addsub = data{addsubline,userInputcol}{:};
if contains(addsub,'y')
    subprefix = 'sub-';
else
    subprefix = '';
end

addsesline = find(strcmp(data.Properties.RowNames,'add ses prefix')); % if session code does not begin with 'ses-' this needs to be prepended for BIDS compatibility
addses = data{addsesline,userInputcol}{:};
if contains(addses,'y')
    sesprefix = 'ses-';
else
    sesprefix = '';
end

% Number of sessions and IDs
dum = find(contains(data.Properties.RowNames,'session ID'));
for i = 1:length(dum)
    ses_id{i} = data{dum(i),userInputcol}{:};
    suff{i} = char(data{dum(i)+1,userInputcol}{:}); % suffix that belongs to this session ID is expected in cell below ID! Use char that in case suffix-cell is empty the routine works with an empthy char
end

% Number of runs/MR series per session; note that program will check in each session. 
% Number of images in each scan series CAUTION: intervals of min-max scans per series must not overlap.
dum = find(contains(data.Properties.RowNames,'MRI series'));
for i = 1:length(dum)
    tasks(i) = data{dum(i),userInputcol};
    taskruns(i) = data{dum(i),runNumbercol};
    metadata_series{i} = char(data{dum(i)+1,userInputcol}{:}); % meta-data info. Use char that in case suffix-cell is empty the routine works with an empty char
    n_series(i,:) = data{dum(i),[minImagescol,maxImagescol]};
    
    if contains(tasks(i),'fieldmap') && contains(tasks(i),'phasediff') 
        
        % Look out for TaskName of corresponding EPI-scan in []
        getTaskName = findstr(tasks{i},'[');        
        if ~isempty(getTaskName)
            getStop = findstr(tasks{i},']');
            
            % Sanity check
            if length(getTaskName)~=length(getStop)
                fprintf('Fieldmap phasediff: wrong usage of square brackets in datasheet.\nProgram stops.\n');
                return
            end
            
            for j = 1:length(getTaskName)
                fmapIntendedFor_series{i,j} = tasks{i}(getTaskName(j)+1:getStop(j)-1);
            end
        end
        
        % B0FieldID in ()
        getFieldID = findstr(tasks{i},'(');        
        if ~isempty(getFieldID)
            getStop = findstr(tasks{i},')');
            
            % Sanity check
            if length(getFieldID)~=length(getStop)
                fprintf('Fieldmap phasediff: wrong usage of brackets in datasheet.\nProgram stops.\n');
                return
            end

            fmapB0FieldID_series{i} = tasks{i}(getFieldID+1:getStop-1);

        else
            fprintf('Fieldmap phasediff: definition of B0FieldIdentifier lacking in datasheet.\n');
            fmapB0FieldID_series{i} = '';

        end


    elseif contains(tasks(i),'dwi') % dwi implemented in v0.4

        % DWI parameter acq is expected in [edgy brackets], DWI parameter dir is expected in (round brackets)
        getAcqName = findstr(tasks{i},'[');        
        if ~isempty(getAcqName)
            getStop = findstr(tasks{i},']');
            
            % Sanity check
            if length(getAcqName)~=length(getStop)
                fprintf('DWI: wrong usage of square brackets in datasheet.\nProgram stops.\n');
                return
            end

                dwi_acq{i} = tasks{i}(getAcqName+1:getStop-1);
        else
            fprintf('DWI: definition of acq parameter lacking.\nProgram stops.\n');
            return
        end

        getDirName = findstr(tasks{i},'(');        
        if ~isempty(getDirName)
            getStop = findstr(tasks{i},')');
            
            % Sanity check
            if length(getDirName)~=length(getStop)
                fprintf('DWI: wrong usage of round brackets in datasheet.\nProgram stops.\n');
                return
            end

                dwi_dir{i} = tasks{i}(getDirName+1:getStop-1);
        else
            fprintf('DWI: definition of acq parameter lacking.\nProgram stops.\n');
            return
        end

    elseif contains(tasks(i),'functional')

        % Task name expected in [edgy brackets], B0 Field Source for B0 Field correction in (round brackets)
        getTaskName = findstr(tasks{i},'[');        
        if ~isempty(getTaskName)
            getStop = findstr(tasks{i},']');
            
            % Sanity check
            if length(getTaskName)~=length(getStop)
                fprintf('Functional: wrong usage of square brackets in datasheet.\nProgram stops.\n');
                return
            end

                func_name{i} = tasks{i}(getTaskName+1:getStop-1);
        else
            fprintf('Functional: definition of task name lacking.\nProgram stops.\n');
            return
        end

        getB0Name = findstr(tasks{i},'(');        
        if ~isempty(getB0Name)
            getStop = findstr(tasks{i},')');
            
            % Sanity check
            if length(getB0Name)~=length(getStop)
                fprintf('Functional: wrong usage of round brackets in datasheet.\nProgram stops.\n');
                return
            end

                B0_sourceID{i} = tasks{i}(getB0Name+1:getStop-1);
        else
            fprintf('Functional: definition of B0FieldSource lacking in datasheet.\n');
            B0_sourceID{i} = '';
        end
    end
end

% Overwrite existing session directories in dataset?
Overwriteline = find(strcmp(data.Properties.RowNames,'overwrite import'));
overwrite = data{Overwriteline,userInputcol}{:};

% Write scan protocol to dicom directory?
WriteProtline = find(strcmp(data.Properties.RowNames,'write scan protocol'));
WriteProt = data{WriteProtline,userInputcol}{:};

logdir = [data_analysis_path,filesep,'RaBIDS-logs'];
if ~isfolder(logdir)
    mkdir(logdir);
end

%% Import data to BIDS structure

% check for generic files and if not yet existant, copy-paste templates
if ~exist(fullfile(data_analysis_path,data_dir,'README'))
    copyfile(fullfile(data_analysis_path,'dataset','code','RaBIDS templates','README'),fullfile(data_analysis_path,data_dir,'README'));
end

if ~exist(fullfile(data_analysis_path,data_dir,'dataset_description.json'))
    copyfile(fullfile(data_analysis_path,'dataset','code','RaBIDS templates','dataset_description.json'),fullfile(data_analysis_path,data_dir,'dataset_description.json'));
end

if ~exist(fullfile(data_analysis_path,data_dir,'.bidsignore'))
    copyfile(fullfile(data_analysis_path,'dataset','code','RaBIDS templates','.bidsignore'),fullfile(data_analysis_path,data_dir,'.bidsignore'));
end

% Import subject data
for i = 1:length(subj_list)
    
    fprintf(['\n\n--------------- Dicom import subject: ',subj_list(i).name{1},' ---------------\n']);
    
    for j = 1:length(ses_id)
        x = clock;
        xs = sprintf('%d-%d-%d-%d-%d-%2.0f',x);
        
        if strcmp(ses_id,'none')
            subject_dir = [data_analysis_path,filesep,data_dir,filesep,subprefix,subj_list(i).name{1}];
            diaryname = [subprefix,subj_list(i).name{1},'_',xs,'_RaBIDS-Import-log.txt'];
        else
            subject_dir = [data_analysis_path,filesep,data_dir,filesep,subprefix,subj_list(i).name{1},filesep,ses_id{j}];
            diaryname = [subprefix,subj_list(i).name{1},'_',sesprefix,ses_id{j},'_',xs,'_RaBIDS-Import-log.txt'];
        end
        
        diary(fullfile(logdir,diaryname));
        
        if ~isfolder(subject_dir) || contains(overwrite,'y')
            fprintf(['\n--------------- Session-ID ',ses_id{j},' ---------------\n\n']);
            
            [out,scanprotocol] = obtain_scanprotocol(HowExpectDicoms,dicomdir,subj_list(i).name{1},suff{j},ses_id{j},general_suffix,data_analysis_path,max_series,addsub,addses,WriteProt);
            [lines,~] = size(out);
            stop = false;
            for nmessages = 1:lines
                fprintf(out{nmessages})
                if contains(out{nmessages},'Program stops')
                    stop = true;
                    break
                end
            end
            if stop
                break
            end

            for k = 1:length(tasks)
                task = tasks{k};
                fprintf(['--------------- Task: ',task,' ---------------\n\n']);
                n_min = n_series(k,1);
                n_max = n_series(k,2);
                metadata = metadata_series{k};
                if contains(task,'fieldmap phasediff')
                    fmapIntendedFor = fmapIntendedFor_series(k,:);
                    B0FieldID = fmapB0FieldID_series{k};
                else
                    fmapIntendedFor = [];
                    B0FieldID = [];
                end
                
                % select MR series from scanprotocol
                [nr,~] = size(scanprotocol);
                valid_series = 'false';
                for l = 1:nr
                    if strcmp(scanprotocol.name{l},metadata)
                        if (scanprotocol.vols(l)>=n_min && scanprotocol.vols(l)<=n_max) || (scanprotocol.vols(l)==n_min && isnan(n_max))
                            useseries = scanprotocol.series(l);
                            valid_series = 'true';
                            break
                        end
                    end
                end
                
                if strcmp(valid_series,'false')
                    fprintf(['No series for task ',task,' with number of volumes MinImages ',num2str(n_min),'-MaxImages ',num2str(n_max),' found.\nWarning #8.\n\n'])
                else                    
                    if strcmp(task,'anat')
                        out = import2nifti_anatomical(HowExpectDicoms,dicomdir,subj_list(i).name{1},suff{j},ses_id{j},general_suffix,data_analysis_path,useseries,scanprotocol.name{useseries},addsub,addses,overwrite,dcm2niix_exe_path);
                    elseif contains(task,'fieldmap')
                        out = import2nifti_fieldmap(HowExpectDicoms,dicomdir,subj_list(i).name{1},suff{j},ses_id{j},general_suffix,data_analysis_path,useseries,scanprotocol.name{useseries},fmapIntendedFor,B0FieldID,addsub,addses,overwrite,dcm2niix_exe_path);
                    elseif contains(task,'dwi')
                        fprintf(['--------------- Run: ',sprintf('%d',taskruns(k)),' ---------------\n\n']);
                        out = import2nifti_diffusion(HowExpectDicoms,dicomdir,subj_list(i).name{1},suff{j},ses_id{j},general_suffix,data_analysis_path,useseries,scanprotocol.name{useseries},dwi_acq{k},dwi_dir{k},taskruns(k),addsub,addses,overwrite,dcm2niix_exe_path);
                    elseif contains(task,'functional')      
                        fprintf(['--------------- Run: ',sprintf('%d',taskruns(k)),' ---------------\n\n']);
                        out = import2nifti_functional(HowExpectDicoms,dicomdir,subj_list(i).name{1},suff{j},ses_id{j},general_suffix,data_analysis_path,useseries,scanprotocol.name{useseries},func_name{k},taskruns(k),B0_sourceID{k},addsub,addses,overwrite,dcm2niix_exe_path);
                    end
                    [lines,~] = size(out);
                    stop = false;
                    for nmessages = 1:lines
                        fprintf(out{nmessages})
                        if contains(out{nmessages},'Program stops')
                            stop = true;
                            break
                        end
                    end
                    if stop
                        break
                    end
                end
            end
        else
            fprintf(['Existing directory found, skip session ID ',ses_id{j},'.\n']);
        end
        diary off
        
    end
end

fprintf('RaBIDS arrived - ready to jump on!\n\n');