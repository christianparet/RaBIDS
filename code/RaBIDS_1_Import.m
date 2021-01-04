%% RaBIDS - Rapid preparation of BIDS
% v0.2.2 release

% Change log 

% v0.2.2 
% - added fieldmap import

% v0.1: 
% - see obtain_scanprotocol; scan protocol saved to subject directory
% - save command window output to RaBIDS-logs directory
% - new entries added to bidsignore-file to allow above output in BIDS validation
% - Datasheet variable name naming convention changed
% - Implemented fieldmap import
% - Error reporting with reference to RaBIDS Error Reference-sheet
% - Removed several bugs

% Change log to v0.2.1
% - for fieldmap phasedifference maps: read TEs from datasheet and push them to fieldmap-nifti-import function

%% Import imaging data to BIDS format
% Program automatically imports your MRI data to nifti and rearranges it according to the BIDS data format

% Christian Paret, ZI Mannheim, 2019-2020; christian.paret@zi-mannheim.de

clear
clc
fprintf('------ RaBIDS ---------- Rapid preparation of BIDS ------- is running ------\n\n');
fprintf('For automated production of .json files during dicom import: open dicm2nii GUI and select option ''save json file''\n');

%% Read data from datasheet
data = readtable('datasheet.xlsx','ReadRowNames',true,'PreserveVariableNames',true,'NumHeaderLines',0);
data_dir =  'dataset';

max_series = 30; % N needs to be equal or larger maximum number of scan series in a subject; this could be improved in future versions

descriptcol = find(strcmp(data.Properties.VariableNames,'Description'));
userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));
minImagescol = find(strcmp(data.Properties.VariableNames,'MinImages'));
maxImagescol = find(strcmp(data.Properties.VariableNames,'MaxImages'));

AddPathline = find(contains(data.Properties.RowNames,'add path'));
for i = 1:length(AddPathline)
    addpath(data{AddPathline(i),userInputcol}{:});
end

DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};

Dicomsline = find(strcmp(data.Properties.RowNames,'dicoms'));
HowExpectDicoms = data{Dicomsline,userInputcol}{:};
switch HowExpectDicoms
    case 'BIDS'
        ExchangeServerPathline = find(strcmp(data.Properties.RowNames,'data exchange path'));
        dicomdir = data{ExchangeServerPathline,userInputcol}{:};
    case 'allinone'
        dicomdir =  'dicomdir';
end

dum = find(contains(data.Properties.RowNames,'subject info'));
for i = 1:length(dum)
    subj_list(i).name = data{dum(i),userInputcol};
end

general_suffixline = find(strcmp(data.Properties.RowNames,'general suffix'));
general_suffix = data{general_suffixline,userInputcol}{:}; % study code at scanner site, used to filter dicoms

addsubline = find(strcmp(data.Properties.RowNames,'add prefix')); % if subject code does not begin with 'sub-' this needs to be prepended for BIDS compatibility
addsub = data{addsubline,userInputcol}{:};
if contains(addsub,'y')
    prefix = 'sub-';
else
    prefix = '';
end

first_imageline = find(strcmp(data.Properties.RowNames,'first image'));
if contains(data{first_imageline,userInputcol}{:},'y')
    first_image = data{first_imageline,minImagescol}; % X=1:[first_image-1] images are skipped and will be deleted from the dicom-directory
else 
    first_image = 0;
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
    metadata_series{i} = char(data{dum(i)+1,userInputcol}{:}); % meta-data info. Use char that in case suffix-cell is empty the routine works with an empty char
    n_series(i,:) = data{dum(i),[minImagescol,maxImagescol]};
    
    if contains(tasks(i),'fieldmap') && contains(tasks(i),'phasediff') 
        
        % Look out for TaskName of corresponding EPI-scan
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
        
        % Read TE1 and TE2 from colomn MinImages and MaxImages
        fmapTE_series(i,:) = data{dum(i)+1,[minImagescol,maxImagescol]};
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
            subject_dir = [data_analysis_path,filesep,data_dir,filesep,prefix,subj_list(i).name{1}];
            diaryname = [prefix,subj_list(i).name{1},'_',xs,'_RaBIDS-Import-log.txt'];
        else
            subject_dir = [data_analysis_path,filesep,data_dir,filesep,prefix,subj_list(i).name{1},filesep,ses_id{j}];
            diaryname = [prefix,subj_list(i).name{1},'_',ses_id{j},'_',xs,'_RaBIDS-Import-log.txt'];
        end
        
        diary(fullfile(logdir,diaryname));
        
        if ~isfolder(subject_dir) || contains(overwrite,'y')
            fprintf(['--------------- Session ID ',ses_id{j},'---------------\n\n']);
            
            [out,scanprotocol] = obtain_scanprotocol(HowExpectDicoms,dicomdir,subj_list(i).name{1},suff{j},ses_id{j},general_suffix,data_analysis_path,max_series,addsub,WriteProt);
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
                if contains(task,'fieldmap')
                    fmapTE1 = fmapTE_series(k,1);
                    fmapTE2 = fmapTE_series(k,2);
                    fmapIntendedFor = fmapIntendedFor_series(k,:);
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
                    fprintf(['No series for task ',task,' with number of volumes MinImages ',num2str(n_min),'-MaxImages ',num2str(n_max),' found.\nError #8.\n\n'])
                else
                    if strcmp(task,'anat')
                        out = import2nifti_anatomical(HowExpectDicoms,dicomdir,subj_list(i).name{1},suff{j},ses_id{j},general_suffix,data_analysis_path,useseries,scanprotocol.name{useseries},addsub,overwrite);
                    elseif contains(task,'fieldmap')
                        out = import2nifti_fieldmap(HowExpectDicoms,dicomdir,subj_list(i).name{1},suff{j},ses_id{j},general_suffix,data_analysis_path,useseries,scanprotocol.name{useseries},fmapTE1,fmapTE2,fmapIntendedFor,addsub,overwrite);
                    else          
                        out = import2nifti_functional(HowExpectDicoms,dicomdir,subj_list(i).name{1},suff{j},ses_id{j},general_suffix,first_image,data_analysis_path,useseries,scanprotocol.name{useseries},task,addsub,overwrite);
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