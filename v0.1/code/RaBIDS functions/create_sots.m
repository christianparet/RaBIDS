function out = create_sots(subject,task,data_analysis_path,ses_id,addsub,firstpulse,condfile,overwrite)
% Create_SOTS_NFB by C. Paret, ZI-Mannheim, 2019-2020
% v0.1 release

%% Comment out if function in use
% clear
% clc
% condfile = 'S:\AG-Austausch\RaBIDS\example\your project directory\dataset\code\conditions_scenes.xlsx';
% data_analysis_path = 'S:\AG-Austausch\RaBIDS\example\your project directory';
% subject = 'RaBIDS01';
% ses_id = 'ses-01';
% addsub = 'y'; % add prefix 'sub-' before subject name
% task = 'scenes';
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
    conddata = readtable(condfile,'ReadRowNames',true);
    
    Namecol = find(strcmp(conddata.Properties.VariableNames,'Name'));
    OnsetIDcol = find(strcmp(conddata.Properties.VariableNames,'OnsetID'));
    Durationcol = find(strcmp(conddata.Properties.VariableNames,'Duration'));
    OffsetIDcol = find(strcmp(conddata.Properties.VariableNames,'OffsetID'));
    
    condlines = find(contains(conddata.Properties.RowNames,'Condition'));
    
    for i = 1:length(condlines)
        allcond{i}.Name = conddata{condlines(i),Namecol}{:};
        allcond{i}.OnsetID = conddata{condlines(i),OnsetIDcol}{:};
        allcond{i}.Duration = conddata{condlines(i),Durationcol}{:};
        allcond{i}.OffsetID = conddata{condlines(i),OffsetIDcol}{:};
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
        logdir = [data_analysis_path,filesep,'sourcedata',filesep,prefix,subject];
        savedir = [data_analysis_path,filesep,'dataset',filesep,prefix,subject,filesep,'func'];
    else
        write_ses = ['_',ses_id,'_'];
        logdir = [data_analysis_path,filesep,'sourcedata',filesep,prefix,subject,filesep,ses_id];
        savedir = [data_analysis_path,filesep,'dataset',filesep,prefix,subject,filesep,ses_id,filesep,'func'];
    end
    
    sotsf = [savedir,filesep,prefix,subject,write_ses,'task-',task,'_multicond.mat'];
    eventsf = [savedir,filesep,prefix,subject,write_ses,'task-',task,'_events.tsv'];
    
    if ~isfile(sotsf) || contains(overwrite,'y') % check whether file exists or overwrite
        
        mkdir(savedir)
        
        %% Read Presentation-logfile, get onsets and durations
        LogFormatLine = find(strcmp(conddata.Properties.RowNames,'Logfile ID format'));
        LogFormat = conddata{LogFormatLine,Namecol}{:};
        LogIDline = find(strcmp(conddata.Properties.RowNames,'Log ID'));
        LogID = conddata{LogIDline,Namecol}{:};
        
        if strcmp(LogFormat,'BIDS')
            logfile = dir([logdir,filesep,prefix,subject,write_ses,LogID,'.log']);
        elseif strcmp(LogFormat,'free')
            logfile = dir([logdir,filesep,LogID,'.log']);
        else
            out{1} = 'User input to object type ''Logfile ID format'' in conditions_TaskName.xlsx is not valid. Select either ''BIDS'' or ''free''.\nProgram stops.\\n';
            return
        end
        
        try
            [event_type,event_code,event_time]=textread([logfile.folder,filesep,logfile.name],'%*d %s %s %d %*s%*s%*s%*s%*s%*s%*s%*s%*s','headerlines',IgnoreHeaderLines);
        catch
            out{1} = ['Cannot read logfile. Check logfile ID definitions in conditions_TasName.xlsx. Check logfile headerlines: data expected to start not before row ',num2str(IgnoreHeaderLines),'.\nContinue with next session.\n\n'];
            return
        end
        
        pulse = 0;
        get_offset = 0;
        
        relevant_types{1}='Pulse';
        relevant_types{2}='Picture';
        relevant_types{3}='Manual';
        
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
                                eval([allcond{i}.Name,'_dur = [',allcond{i}.Name,'_dur str2num(allcond{i}.Duration) ];']); % ... append duration vector
                                get_offset = 0;
                            else
                                get_offset = 1; % If duration is not defined (i.e. cell in conditions-table empty)...
                                
                                if allcond{i}.OffsetID >= 0 % ... and offset ID is defined ...
                                    retain_offsetCondName = allcond{i}.Name; % ... then keep some information in memory to define duration dependent from offset ID
                                    retain_offsetID = allcond{i}.OffsetID;
                                    on_time = (event_time(k)-t0)/10000;
                                    out{dum,:} = ['No duration defined for condition: ', allcond{i}.Name,'. OffsetID is used.\n'];
                                    dum = dum + 1;
                                else % Duration or offset ID must be defined, otherwise: force function to stop
                                    out{dum,:} = ['Duration and OffsetID definition missing for condition ', allcond{i}.Name,'.\nProgram stops.\n'];
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
        
        for i = 1:length(allcond)
            
            % spm-readable format
            names{i} = allcond{i}.Name;
            eval(['onsets{i} = ',allcond{i}.Name,'_onset ;']);
            eval(['durations{i} = ',allcond{i}.Name,'_dur ;']);
            
           % BIDS format
            all_trialtypes = horzcat(all_trialtypes, repmat({allcond{i}.Name},[1 length(onsets{i})]));
            eval(['all_onsets = [all_onsets,',allcond{i}.Name,'_onset] ;']);
            eval(['all_durations = [all_durations,',allcond{i}.Name,'_dur] ;']);
            
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
        out{dum,:} = 'SOTs files exist, continue with next step.\n';
    end
    
catch
    out{dum,:} = 'An error ouccured. Program stops\n';
end
