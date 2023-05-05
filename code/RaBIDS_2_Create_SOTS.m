%% RaBIDS - Rapid analysis with BIDS
% Most recently updated for v0.3 release

%% Create SOTS
% Chris Paret, ZI Mannheim, 2019-2021; christian.paret@zi-mannheim.de
% This program preprocesses and (first-level-) analyzes fMRI data
% Run RaBIDS_Import first to create BIDS data structure
% Edit datasheet and conditions_TASK.xlsx file (TASK = your task name); see explanation in file

% Change log from v0.1/v0.2:
% - bugfix to enable create sots with the delete-option of initial MRI images being switched off
% - save command window output to logfile directory
% - datasheet variable naming convention changed
% - sourcedata path from datasheet

clear
clc
fprintf('------ RaBIDS ---------- Rapid analysis with BIDS ------- is running ------\n\n');

%% Define task and text-file with definitions
task = input('Type ''TaskName''\n'); % enter task name as defined in datasheet
condfile = ['conditions_',task,'.xlsx']; 

%% Read data from datasheet
% Entries of this paragraph need to be edited by user for fMRI-task to analyze:
data = readtable('datasheet.xlsx','ReadRowNames',true,'PreserveVariableNames',true,'NumHeaderLines',0);
data_dir =  'dataset';

descriptcol = find(strcmp(data.Properties.VariableNames,'Description'));
userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));
minImagescol = find(strcmp(data.Properties.VariableNames,'MinImages'));

AddPathline = find(contains(data.Properties.RowNames,'add path'));
for i = 1:length(AddPathline)
    addpath(data{AddPathline(i),userInputcol}{:});
end

DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};

SourceDataPathline = find(strcmp(data.Properties.RowNames,'sourcedata path'));
source_data_path = data{SourceDataPathline,userInputcol}{:};

dum = find(contains(data.Properties.RowNames,'subject info'));
for i = 1:length(dum)
    subj_list(i).name = data{dum(i),userInputcol};
end

addsubline = find(strcmp(data.Properties.RowNames,'add prefix')); % if subject code does not begin with 'sub-' this needs to be prepended for BIDS compatibility
addsub = data{addsubline,userInputcol}{:};
if contains(addsub,'y')
    prefix = 'sub-';
else
    prefix = '';
end

first_imageline = find(strcmp(data.Properties.RowNames,'first image'));
first_image = data{first_imageline,minImagescol}; % this item identifies timepoint t0 based on pulse number
if first_image > 1
    fprintf(['Stimulus timings are defined relative to volume nr ', num2str(first_image),'.\n']);
else 
    first_image = 1;
    fprintf('Stimulus timings are defined relative to the first volume available.\n');
end

if contains(data{first_imageline,userInputcol}{:},'y')
    fprintf(['According to user input MRI images up to volume nr ',num2str(first_image),' have been deleted.\nStimulus timings should be in sync with your volumes.\n']);
else
    fprintf(['According to user input initial MRI images were not removed.\nIt will be critical to discard images up to volume nr ',num2str(first_image),' in firstlevel GLM analysis.\n']);
end

% Number of sessions and IDs
dum = find(contains(data.Properties.RowNames,'session ID'));
for i = 1:length(dum)
    ses_id{i} = data{dum(i),userInputcol}{:};
    suff{i} = char(data{dum(i)+1,userInputcol}{:}); % suffix that belongs to this session ID is expected in cell below ID! Use char that in case suffix-cell is empty the routine works with an empthy char
end

% Overwrite existing session directories in dataset
Overwriteline = find(strcmp(data.Properties.RowNames,'overwrite sots'));
overwrite = data{Overwriteline,userInputcol}{:};

% Make directory to save RaBIDS logfiles
logdir = [data_analysis_path,filesep,'RaBIDS-logs'];
if ~isfolder(logdir)
    mkdir(logdir);
end

%% Create sots
% Developed for Presentation logfiles (Neurobehavioral Systems)
fprintf(['\n------ Create SOTs -------\nTask: ',task,'\nFile: ',condfile]);
for i = 1:length(subj_list)
    fprintf(['\nAnalyze subject: ',subj_list(i).name{1},'\n']);
    
    for j = 1:length(ses_id)
        x = clock;
        xs = sprintf('%d-%d-%d-%d-%d-%2.0f',x);
        
        if strcmp(ses_id,'none')
            subject_dir = [data_analysis_path,filesep,data_dir,filesep,prefix,subj_list(i).name{1}];
            diaryname = [prefix,subj_list(i).name{1},'_',xs,'_RaBIDS-Create_SOTS-log.txt'];
        else
            subject_dir = [data_analysis_path,filesep,data_dir,filesep,prefix,subj_list(i).name{1},filesep,ses_id{j}];
            diaryname = [prefix,subj_list(i).name{1},'_',ses_id{j},'_',xs,'_RaBIDS-Create_SOTS-log.txt'];
        end
        
        diary(fullfile(logdir,diaryname));
        
        fprintf(['\nSession ID: ',ses_id{j},'\n']);
        out = create_sots(subj_list(i).name{1},task,data_analysis_path,source_data_path,ses_id{j},addsub,first_image,condfile,overwrite);
        [lines,~] = size(out);
        stop = false;
        for k = 1:lines
            fprintf(out{k})
            if contains(out{k},'Program stops\n')
                stop = true;
                break
            end
        end
        diary off
        
        if stop
            break
        end
    end
end

fprintf('RaBIDS arrived - ready to jump on!\n\n');