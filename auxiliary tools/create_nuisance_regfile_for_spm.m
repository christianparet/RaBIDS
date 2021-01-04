%% Creates nuisance regressor mat-file ready for spm analysis
% Christian Paret, ZI-Mannheim, 2021

% Requires dataset in BIDS format and derivatives output from fmriPrep (tested with version 20.2.0)
% 6 realignment regressors (3 translation, 3 rotation parameters) are used
% Download program to dataset/code directory to run

%% Read data from datasheet
data = readtable('datasheet.xlsx','ReadRowNames',true,'PreserveVariableNames',true,'NumHeaderLines',0);

userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));

DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};

%% identify derivative directory to work on
derivd = dir([data_analysis_path,filesep,'derivatives',filesep,'fmriprep*']);

fprintf(['Found ',num2str(length(derivd)),' derivative directories:\n'])

for i=1:length(derivd)
    fprintf([num2str(i),') ',derivd(i).name,'\n'])
end

selectd = input('Enter number of derivative directory to unzip preporcessed niftis and press enter.\n');
derivd = [data_analysis_path,filesep,'derivatives',filesep,derivd(selectd).name];

%% Task to work on
reqtask = input('What task to work on? Enter a single task name or enter ''all'' to process all tasks available.\n');

%% Load confound timeseries produced by fMRIPrep and save spm nuisance regressor files
% in case of nested structure: account for 2 layers of fmriprep directories
try
    childd = dir([derivd,filesep,'fmriprep*']);
    if length(childd)>1
        fprintf('Too many subdirectories.\n')
        return
    else
        derivd = [derivd,filesep,childd.name];
    end
end

try
    childd = dir([derivd,filesep,'fmriprep*']);
    if length(childd)>1
        fprintf('Too many subdirectories.\n')
        return
    else
        derivd = [derivd,filesep,childd.name];
    end
end

%% Read confounds timeseries file and write realignment parameters to nuisance file

allsubs = dir([derivd,filesep,'sub-*']);

for sub = 1:length(allsubs)
    subd = [derivd,filesep,allsubs(sub).name];
    
    if isfolder(subd)
        
        fprintf([allsubs(sub).name,'\n'])
        
        try
            allses = dir([subd,filesep,'ses-*']);
            
            for ses = 1:length(allses)
                sesd = [subd,filesep,allses(ses).name];
                fprintf([allses(ses).name,'\n'])
                
                if isfolder([sesd,filesep,'func'])
                    conf_timeseriesf = dir([sesd,filesep,'func',filesep,'sub-*-confounds_timeseries.tsv']);
                    
                    for funcs = 1:length(conf_timeseriesf)
                        pos1 = strfind(conf_timeseriesf(funcs).name,'task');
                        pos2 = strfind(conf_timeseriesf(funcs).name,'_desc-confounds_timeseries');
                        taskid = conf_timeseriesf(funcs).name(pos1+5:pos2-1);
                        
                        if strcmp(reqtask,'all') || strcmp(reqtask,taskid)
                            conf_table = readtable([sesd,filesep,'func',filesep,conf_timeseriesf(funcs).name],'FileType','text','ReadRowNames',false,'PreserveVariableNames',true,'NumHeaderLines',0);
                            R = [conf_table.trans_x, conf_table.trans_y, conf_table.trans_z, conf_table.rot_x, conf_table.rot_y, conf_table.rot_z];
                            save(fullfile(sesd,'func',[allsubs(sub).name,'_',allses(ses).name,'_task-',taskid,'_desc-confounds_timeseries_desc-realignment_regressors.mat']),'R')
                            fprintf(['Confounds found for ',conf_timeseriesf(funcs).name,'.\n'])
                        end
                    end
                else
                    fprintf('No func directory found or confound file not readable.\n');
                end
            end
        catch
            if isfolder([subd,filesep,'func'])
                conf_timeseriesf = dir([subd,filesep,'func',filesep,'sub-*-confounds_timeseries.tsv']);
                
                for funcs = 1:length(conf_timeseriesf)
                    pos1 = strfind(conf_timeseriesf(funcs).name,'task');
                    pos2 = strfind(conf_timeseriesf(funcs).name,'_desc-confounds_timeseries');
                    taskid = conf_timeseriesf(funcs).name(pos1+5:pos2-1);

                    if strcmp(reqtask,'all') || strcmp(reqtask,taskid)
                        conf_table = readtable([subd,filesep,'func',filesep,conf_timeseriesf(funcs).name],'FileType','text','ReadRowNames',false,'PreserveVariableNames',true,'NumHeaderLines',0);
                        R = [conf_table.trans_x, conf_table.trans_y, conf_table.trans_z, conf_table.rot_x, conf_table.rot_y, conf_table.rot_z];
                        save(fullfile(subd,'func',[allsubs(sub).name,'_task-',reqtask,'_desc-confounds_timeseries_desc-realignment_regressors.mat']),'R')
                        fprintf(['Confounds found for ',conf_timeseriesf(funcs).name,'.\n'])
                    end
                end
            else
                fprintf('No func directory found.\n');
            end
        end
    end
end
fprintf('End\n')