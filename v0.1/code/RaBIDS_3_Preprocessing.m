%% RaBIDS - Rapid analysis with BIDS
% v0.1 release

%% Preprocessing step
% Chris Paret, ZI Mannheim, 2019-2020; christian.paret@zi-mannheim.de
% This program preprocesses fMRI data
% Run RaBIDS_Import and RaBIDS_Create_SOTS first to create BIDS data structure
% Edit datasheet and conditions_TASK.xlsx file (TASK = your task name); see explanation in file

clear
clc
fprintf('------ RaBIDS ---------- Rapid analysis with BIDS ------- is running ------\n\n');

%% Define task and text-file with definitions
task = input('Type ''TaskName''\n'); % enter task name as defined in datasheet

%% Read data from datasheet
% Entries of this paragraph need to be edited by user for fMRI-task to analyze:
data = readtable('datasheet.xlsx','ReadRowNames',true);

descriptcol = find(strcmp(data.Properties.VariableNames,'description'));
userInputcol = find(strcmp(data.Properties.VariableNames,'userInput'));
minImagescol = find(strcmp(data.Properties.VariableNames,'minimumImages'));

AddPathline = find(contains(data.Properties.RowNames,'add path'));
for i = 1:length(AddPathline)
    addpath(data{AddPathline(i),userInputcol}{:});
end

DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};

dum = find(contains(data.Properties.RowNames,'subject info'));
for i = 1:length(dum)
    subj_list(i).name = data{dum(i),userInputcol};
end

addsubline = find(strcmp(data.Properties.RowNames,'add prefix')); % if subject code does not begin with 'sub-' this needs to be prepended for BIDS compatibility
addsub = data{addsubline,userInputcol}{:};

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

% Overwrite existing session directories in dataset
Overwriteline = find(strcmp(data.Properties.RowNames,'overwrite preprocess'));
overwrite = data{Overwriteline,userInputcol}{:};

%% Preprocess data
% Preprocessed data are written to derivatives directory
fprintf(['\n------ Preprocess fMRI data -------\nTask: ',task]);
[spmpath,~,~] = fileparts(which('spm'));
for i = 1:length(subj_list)
    fprintf(['\nAnalyze subject: ',subj_list(i).name{1},'\n']);
    
    for j = 1:length(ses_id)
        fprintf(['Session ID: ',ses_id{j},'\n']);
        out = prepro(subj_list(i).name{1},task,data_analysis_path,spmpath,ses_id{j},addsub,overwrite);
        fprintf(out)
        if contains(out,'Program stops\n')
            break
        end
    end
end

fprintf('RaBIDS arrived - ready to jump on!\n\n');