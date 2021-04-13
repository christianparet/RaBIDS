%% Find spikes and trim session
% Christian Paret, ZI Mannheim, 2021
%
% CAVE: current version is developed for efp-study and may not generalize to other studies
%
% The script assists in detecting outlier motions and trimming sessions before/after outlier motions. This can be useful if large motions are detected in the beginning or the end of a session, and if removal of images before/after that timepoint is not expected to lead to major data loss.
% Before working with this script, framewise displacement parameters should be visually inspected for above-threshold motion. Affected sessions and tasks should be known by the user for efficient use.

clear
clc

%% User input required
fd_thresh = 4; % threshold for motion outlier based on framewise displacement parameter

%% Read data from datasheet
data = readtable('datasheet.xlsx','ReadRowNames',true,'PreserveVariableNames',true,'NumHeaderLines',0);

userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));

DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};

datasetd = [data_analysis_path,filesep,'dataset'];

%% Identify derivative directory to work on
derivd = dir([data_analysis_path,filesep,'derivatives',filesep,'fmriprep*']);

fprintf(['Found ',num2str(length(derivd)),' derivative directories:\n'])

for i=1:length(derivd)
    fprintf([num2str(i),') ',derivd(i).name,'\n'])
end

selectd = input('Enter number of derivative directory to work on and press enter.\n');
derivd = [data_analysis_path,filesep,'derivatives',filesep,derivd(selectd).name];

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

%% What subject to work on
reqsub = input('What subject to work on?\n');

%% What session to work on
reqses = input('What session to work on?\n');

%% Task to work on
reqtask = input('What task to work on?\n');

%% Find motion outliers
funcd = fullfile(derivd,reqsub,reqses,'func');
conf_timeseriesf = [reqsub,'_',reqses,'_task-',reqtask,'_desc-confounds_timeseries.tsv'];

if isfile([funcd,filesep,conf_timeseriesf])
    conf_table = readtable([funcd,filesep,conf_timeseriesf],'FileType','text','ReadRowNames',false,'PreserveVariableNames',true,'NumHeaderLines',0);
    
    close
    plot(conf_table.framewise_displacement)
    hold on
    plot(1:length(conf_table.framewise_displacement),repmat(fd_thresh,length(conf_table.framewise_displacement))) % plot threshold as line
    title([reqsub,' ',reqses,' ',reqtask]);
    hold off
    
    spikes = find(conf_table.framewise_displacement>fd_thresh);
    if spikes
        fprintf(['Found ',num2str(length(spikes)),' spike(s):\n'])
        for i = 1:length(spikes)
            fprintf([' - Spike nr ',num2str(i,'%02d'),' at ',num2str(spikes(i),'%03d'),'.\n'])
        end
        trim = input('Trim session (Y/N)?\n');
        if strcmp(trim,'Y')
            trimbegin = input('1.) Enter spike nr after which to keep data. If you want to keep all data before a spike enter 0.');
            trimend = input('2.) Enter spike nr after which to discard data. If you want to keep all data after a spike enter 0.');
            
            if ~trimbegin && ~trimend
                fprintf('Your input reads as if you don''t want to discard any data.\n')
                return
            end
            
            if trimbegin
                trimses(1) = spikes(trimbegin)+6;
            else
                trimses(1) = 1;
            end
            
            if trimend
                trimses(2) = spikes(trimend)-6;
            else 
                trimses(2) = length(conf_table.framewise_displacement);
            end
            
            save([funcd,filesep,reqsub,'_',reqses,'_task-',reqtask,'_trimsession'],'trimses')
            fprintf('Note: 5 more scans are discarded around the spike to prevent spurious signals.\nTrimsession file saved to derivatives directory.\n')
        end
    else
        fprintf('No spikes detected.\n')
    end    
        
else
    fprintf('Could not find confounds-timeseries file for this subject/session/task.\n')
end
   
                
