function out = create_sots(subject,task,data_analysis_path,source_data_path,ses_id,addsub,firstpulse,condfile,overwrite)
% Create_SOTS by C. Paret, ZI-Mannheim, 2019-2023
% Most recently updated for v0.3 release

% 2021/08/09: Import options for readtable function revisited to work in a more generic way
% 2022/03 & 2023/04: Import options for readtable function revisited to work in a more generic way: Mirus Jindrová added opts definitions
% 2023/04: Option to split task condition in early and late events. If selected, SOTs-files have two stimulus onset functions.

%% Comment out if function in use
% clear
% clc
% condfile = 'E:\mytrainingdata\your project directory\dataset\code\conditions_faces.xlsx';
% data_analysis_path = 'E:\mytrainingdata\your project directory';
% source_data_path = 'E:\mytrainingdata\data exchange server\RABIDS-example\sourcedata';
% subject = 'RaBIDS01';
% ses_id = 'ses-01';
% addsub = 'y';
% task = 'faces';
% firstpulse = 5;
% overwrite = 'y';
%%
IgnoreHeaderLines = 5; % Number of headerlines to ignore when reading logfile.
dum = 1;

if ~firstpulse
    out{1} = 'User input of object type ''first image'' in ''datasheet.xlsx'' must be ''yes'' and minimum images must be >0';
    return
end

try
    %% Define conditions
    opts = detectImportOptions(condfile,'DataRange','B2','RowNamesRange','A2','VariableNamesRange','B1');
    opts = setvartype(opts, ["Name","OnsetID","OffsetID","ContrastType","ContrastPlus1","ContrastMinus1"],"string");
    opts = setvartype(opts, ["Duration","SplitCondition"],"double");
    conddata = readtable(condfile,opts);
    
    Namecol = find(strcmp(conddata.Properties.VariableNames,'Name'));
    OnsetIDcol = find(strcmp(conddata.Properties.VariableNames,'OnsetID'));
    Durationcol = find(strcmp(conddata.Properties.VariableNames,'Duration'));
    OffsetIDcol = find(strcmp(conddata.Properties.VariableNames,'OffsetID'));
    SplitCondcol = find(strcmp(conddata.Properties.VariableNames,'SplitCondition'));
    
    condlines = find(contains(conddata.Properties.RowNames,'Condition'));
        
    for i = 1:length(condlines)
        allcond{i}.Name = conddata{condlines(i),Namecol}{:};
        allcond{i}.OnsetID = conddata{condlines(i),OnsetIDcol}{:};
        allcond{i}.SplitCond = conddata{condlines(i),SplitCondcol};
                
        if ismissing(conddata{condlines(i),Durationcol})
            if ismissing(conddata{condlines(i),OffsetIDcol})
                out{dum,:} = ['Duration and OffsetID not defined for condition ',allcond{i}.Name,'.\nWarning #4.\nProgram stops\n'];
                return
            else
                allcond{i}.OffsetID = conddata{condlines(i),OffsetIDcol}{:};
                allcond{i}.Duration = -1;
                out{dum,:} = ['Duration not defined for condition ',allcond{i}.Name,'. OffsetID is used.\n'];
                dum=dum+1;
            end
        else
            allcond{i}.Duration = conddata{condlines(i),Durationcol};
        end
        
        eval([allcond{i}.Name,'_onset = [];']); % create empty variable for each condition
        eval([allcond{i}.Name,'_dur = [];']);
    end
    
    %% Define paths
    if contains(addsub,'y')
        prefix = 'sub-';
    else
        prefix = '';
    end
        
    if strcmp(ses_id,'none')
        write_ses = '_';
        logdir = [source_data_path,filesep,prefix,subject];
        
        % give flexibility if files in sourcedata folder do not stringently follow naming scheme
        if ~isfolder(logdir)
            logdir = [source_data_path,filesep,subject];
            if ~isfolder(logdir)
                out{dum,:} = 'Logfile directory not found.\Warning #16\nContinue with next session.\n';
                return
            end
        end
            
        savedir = [data_analysis_path,filesep,'dataset',filesep,prefix,subject,filesep,'func'];
    else
        write_ses = ['_',ses_id,'_'];
        logdir = [source_data_path,filesep,prefix,subject,filesep,ses_id];
        
        % give flexibility if files in sourcedata folder do not stringently follow naming scheme
        if ~isfolder(logdir)
            logdir = [source_data_path,filesep,subject,filesep,ses_id];
            if ~isfolder(logdir)
                out{dum,:} = 'Logfile directory not found.\nWarning #16\nContinue with next session.\n';
                return
            end
        end
        
        savedir = [data_analysis_path,filesep,'dataset',filesep,prefix,subject,filesep,ses_id,filesep,'func'];
    end
    
    sotsf = [savedir,filesep,prefix,subject,write_ses,'task-',task,'_multicond.mat'];
    eventsf = [savedir,filesep,prefix,subject,write_ses,'task-',task,'_events.tsv'];
    
    if ~isfile(sotsf) || contains(overwrite,'y') % check whether file exists or overwrite
        
        if ~isfolder(savedir)
            mkdir(savedir)
        end
        
        %% Read Presentation-logfile, get onsets and durations
        LogFormatLine = find(strcmp(conddata.Properties.RowNames,'Logfile ID format'));
        LogFormat = conddata{LogFormatLine,Namecol}{:};
        LogIDline = find(strcmp(conddata.Properties.RowNames,'Log ID'));
        LogID = conddata{LogIDline,Namecol}{:};
        
        if strcmp(LogFormat,'BIDS')
            filename = [subject,write_ses,LogID,'.log'];
            logfile = dir([logdir,filesep,filename]); % Tolerate deviation of naming scheme           
            if isempty(logfile)
                filename = [prefix,filename];
                logfile = dir([logdir,filesep,filename]);
            else
                out{dum,:} = 'Logfile not found.\nWarning #15\nContinue with next session.\n';
                return
            end
        elseif strcmp(LogFormat,'free')
            filename = [LogID,'.log'];
            logfile = dir([logdir,filesep,filename]);
            if isempty(logfile)
                out{dum,:} = 'Logfile not found.\nWarning #15\nContinue with next session.\n';
                return
            end          
        else
            out{dum,:} = 'User input to object type ''Logfile ID format'' in conditions_TaskName.xlsx is invalid.\nWarning #5.\nProgram stops.\\n';
            return
        end
        
        try
            [event_type,event_code,event_time] = textread([logfile.folder,filesep,logfile.name],'%*d %s %s %d %*s%*s%*s%*s%*s%*s%*s%*s%*s','headerlines',IgnoreHeaderLines);
            out{dum,:} = ['Use logfile ',logfile.name,'\n'];
            dum=dum+1;
        catch
            out{dum,:} = 'Cannot read logfile.\nWarning #6.\nContinue with next session.\n\n';
            return
        end
        
        pulse = 0;
        get_offset = 0;
        
        relevant_types{1}='Pulse';
        relevant_types{2}='Picture';
        relevant_types{3}='Manual';
        relevant_types{4}='Response';
        
        retain_offsetID = 'none';
        
        for k=1:length(event_type)
            
            if any(strcmp(event_type{k},relevant_types))
                
                if strcmp(event_type{k},'Pulse')
                    pulse = pulse + 1;
                    if pulse==firstpulse
                        t0=event_time(k);
                    end
                    
                else
                    for i = 1:length(allcond)
                        if get_offset && contains(event_code{k},retain_offsetID) % This statement is entered if the duration of the previous condition is defined by an OffsetID
                            eval([retain_offsetCondName,'_dur = [',retain_offsetCondName,'_dur (event_time(k)-t0)/10000-on_time ];']);
                            get_offset = 0;
                        end
                        
                        if contains(event_code{k},allcond{i}.OnsetID) % If onset ID from one of the defined conditions is recognized, enter this statement
                            eval([allcond{i}.Name,'_onset = [',allcond{i}.Name,'_onset (event_time(k)-t0)/10000 ];']); % Append time-vector with this event; subtract t0 timestamp
                            
                            if allcond{i}.Duration >= 0 % If duration is defined...
                                eval([allcond{i}.Name,'_dur = [',allcond{i}.Name,'_dur allcond{i}.Duration ];']); % ... append duration vector
                                get_offset = 0;
                            else
                                get_offset = 1; % If duration is not defined (i.e. cell in conditions-table empty)...
                                
                                if ischar(allcond{i}.OffsetID) % ... and offset ID is defined ...
                                    retain_offsetCondName = allcond{i}.Name; % ... then keep some information in memory to define duration dependent from offset ID
                                    retain_offsetID = allcond{i}.OffsetID;
                                    on_time = (event_time(k)-t0)/10000;
%                                     out{dum,:} = ['No duration defined for condition: ', allcond{i}.Name,'. OffsetID is used.\n'];
%                                     dum = dum + 1;
                                else % Duration or offset ID must be defined, otherwise: force function to stop
                                    out{dum,:} = ['Duration and OffsetID definition missing for condition ', allcond{i}.Name,'.\nWarning #4.\nProgram stops.\n'];
                                    dum = dum + 1;
                                    return
                                end
                            end
                        end
                    end
                end
            end
        end
        
        %% Write stimulus onset function
        
        all_trialtypes = [];
        all_onsets = [];
        all_durations = [];

        cond_i = 1;
        
        for i = 1:length(allcond)

            % produce two stimulus functions (i.e., early, late) per condition if required
            if allcond{i}.SplitCond == 2 

                out{dum,:} = ['Condition ',allcond{i}.Name,': SplitCondition selected. RaBIDS produces two stimulus onset functions for early and late events.\n'];
                dum=dum+1;
                
                eval(['binsize = floor ( length ( ',allcond{i}.Name,'_onset ) / 2 ) ; ']);
                
                % spm-readable format
                names{cond_i} = [allcond{i}.Name,'_early'];
                eval(['onsets{cond_i} = ',allcond{i}.Name,'_onset(1:binsize) ;']);
                eval(['durations{cond_i} = ',allcond{i}.Name,'_dur(1:binsize) ;']);
                
                % BIDS format
                all_trialtypes = horzcat(all_trialtypes, repmat({[allcond{i}.Name,'_early']},[1 length(onsets{i})]));
                eval(['all_onsets = [all_onsets,',allcond{i}.Name,'_onset(1:binsize)] ;']);
                eval(['all_durations = [all_durations,',allcond{i}.Name,'_dur(1:binsize)] ;']);

                cond_i = cond_i + 1;

                % spm-readable format
                names{cond_i} = [allcond{i}.Name,'_late'];
                eval(['onsets{cond_i} = ',allcond{i}.Name,'_onset(end-binsize+1:end) ;']);
                eval(['durations{cond_i} = ',allcond{i}.Name,'_dur(end-binsize+1:end) ;']);

                % BIDS format
                all_trialtypes = horzcat(all_trialtypes, repmat({[allcond{i}.Name,'_late']},[1 length(onsets{cond_i})]));
                eval(['all_onsets = [all_onsets,',allcond{i}.Name,'_onset(end-binsize+1:end)] ;']);
                eval(['all_durations = [all_durations,',allcond{i}.Name,'_dur(end-binsize+1:end)] ;']);

                cond_i = cond_i + 1;
            else

                if allcond{i}.SplitCond > 2
                    out{dum,:} = ['Condition ',allcond{i}.Name,': Value set for SplitCondition in conditions_task.xlsx file is > 2, which is invalid. RaBIDS ignores such values.\n'];
                    dum=dum+1;
                end

                % spm-readable format
                names{cond_i} = allcond{i}.Name;
                eval(['onsets{cond_i} = ',allcond{i}.Name,'_onset ;']);
                eval(['durations{cond_i} = ',allcond{i}.Name,'_dur ;']);
            
                % BIDS format
                all_trialtypes = horzcat(all_trialtypes, repmat({allcond{i}.Name},[1 length(onsets{cond_i})]));
                eval(['all_onsets = [all_onsets,',allcond{i}.Name,'_onset] ;']);
                eval(['all_durations = [all_durations,',allcond{i}.Name,'_dur] ;']);

                cond_i = cond_i + 1;
            end
        end
        
        % save spm-readable
        save(sotsf,'names','onsets','durations');
        out{dum,:} = 'SOTs created and saved\n';
        
        % save BIDS readable
        onset = all_onsets';
        duration = all_durations';
        trial_type = all_trialtypes';
        events = table(onset,duration,trial_type);
        writetable(events,eventsf,'Delimiter','tab','FileType','text')        

    else
        out{dum,:} = 'SOTs file exists and permission to overwrite is false.\nContinue with next step.\n';
    end
    
catch
    out{dum,:} = 'Warning #3.\nProgram stops\n';
end
