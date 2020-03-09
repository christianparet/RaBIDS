%% RaBIDS - Rapid analysis with BIDS
% v0.1 release

%% Firstlevel step, uses SPM12
% Chris Paret, ZI Mannheim, 2019-2020; christian.paret@zi-mannheim.de
% This program analyzes fMRI data
% Edit datasheet.xlsx and conditions_TaskName.xlsx file; see explanation in file

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

%% Scrubbing
% Add nuisance regressors to control for super-threshold movements
fprintf(['\n------ Scrubbing -------\nTask: ',task]);
for i = 1:length(subj_list)
    fprintf(['\nAnalyze subject: ',subj_list(i).name{1},'\n']);
    
    for j = 1:length(ses_id)
        fprintf(['Session ID: ',ses_id{j},'\n']);
        out = scrubbing(subj_list(i).name{1},task,data_analysis_path,ses_id{j},addsub,overwrite);
        fprintf(out)
        if contains(out,'Program stops\n')
            break
        end
    end
end

%% First level
fprintf(['\n------ First-level fMRI analysis -------\nTask: ',task]);
condfile = ['conditions_',task,'.xlsx']; 
for i = 1:length(subj_list)
    fprintf(['\nAnalyze subject: ',subj_list(i).name{1},'\n']);
    
    for j = 1:length(ses_id)
        fprintf(['Session ID: ',ses_id{j},'\n']);
        out = firstlevel(subj_list(i).name{1},task,data_analysis_path,ses_id{j},addsub,condfile,overwrite);
        fprintf(out)
        if contains(out,'Program stops\n')
            break
        end
    end
end

fprintf('RaBIDS arrived - ready to jump on!\n\n');