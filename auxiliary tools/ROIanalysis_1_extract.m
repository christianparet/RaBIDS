%% Extract ROI data using the marsbar toolbox
% Christian Paret, ZI Mannheim, 2021

% The script uses the marsbar toolbox to extract parameter estimates, percent signal change and time course data from regions of interest.
% Script can be used on firstlevel files produced with SPManalysis_3_firstlevel.m

% It is ok to ignore matlab Warning message "Action '!ContextHelp' is deprecated.

clear

%% Set path
% Add path to marsbar subfolders
addpath('E:\marsbar-0.44\marsbar-0.44');
folder = fileparts(which('marsbar'));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%% Define directories and initiate variables
pos = 0;
firstleveldir = uigetdir(pwd,'Select firstlevel directory including the subject directories'); % directory containing child directories with each child containing an estimated firstlevel SPM.mat
% firstleveldir = 'Y:\Projects\EFPTest\Data_analysis\spm_analysis\firstlevel\ses-pre\task-efptest\taskrelated_activity';
subdirs = dir([firstleveldir,filesep,'sub*']);
maskdir = uigetdir(pwd,'Select directory containing binary brain mask files (*roi.mat)');
% maskdir = 'Y:\Projects\EFPTest\Data_analysis\brain masks';

%%
if ~isempty(subdirs)
    maskfile = uigetfile(maskdir,'Select a mask file'); % needs to be a mat file
    roif = fullfile(maskdir,maskfile);
    
    dum = strfind(maskfile,'.mat');
    
    %% Define filenames to write results
    filename = fullfile(firstleveldir,['ROI-cons_mask-',maskfile(1:dum),'txt']);
    con_file = fopen(filename,'w');
    filename = fullfile(firstleveldir,['ROI-percent-BOLD-signal-change_mask-',maskfile(1:dum),'txt']);
    psc_file = fopen(filename,'w');
    filename = fullfile(firstleveldir,['ROI-FIR-timecourse_mask-',maskfile(1:dum),'txt']);
    fir_tc_file = fopen(filename,'w');
    filename = fullfile(firstleveldir,['ROI-fitted-timecourse_mask-',maskfile(1:dum),'txt']);
    fitted_tc_file = fopen(filename,'w');
    
    %% Iterate through subjects
    for sub = 1:length(subdirs)
        if subdirs(sub).isdir
            fprintf(['\nProcessing ',subdirs(sub).name,'.\n'])
            dum1 = strfind(subdirs(sub).name,'_ses-');
            if isempty(dum1)
                dum1 = strfind(subdirs(sub).name,'_task-');
                if isempty(dum1)
                    fprintf('An error occured. Directory may not be consistent with BIDS.\n')
                    return
                end
            end
            
            % Duration for which FIRs should be estimated (in seconds)
            fir_length = 24;
            
            % Duration of stimulus events (in seconds)
            e_length = 1;
            
            spm_name = fullfile(firstleveldir,subdirs(sub).name,'SPM.mat');
            
            %% Get contrast values (from marsbar manual: http://marsbar.sourceforge.net/index.html)
            fprintf('Get contrast values.\n')
            % Make marsbar design object
            D = mardo(spm_name);
            % Make marsbar ROI object
            R = maroi(roif);
            % Fetch data into marsbar data object
            Y = get_marsy(R, D, 'mean');
            % Get contrasts from original design
            xCon = get_contrasts(D);
            % Estimate design on ROI data
            E = estimate(D, Y);
            
            % Put contrasts from original design back into design object
            E = set_contrasts(E, xCon);
            % get design betas
            b = betas(E);
            % get stats and stuff for all contrasts into statistics structure
            marsS = compute_contrasts(E, 1:length(xCon));
            
            %% Extract percent signal change (from marsbar manual)
            fprintf('Extract percent signal change.\n')
            %% What "events" are returned below???
            [e_specs, e_names] = event_specs(D);
            n_events = size(e_specs, 2);
            dur = 0;
            % Return percent signal esimate for all events in design
            for e_s = 1:n_events
                pct_ev(e_s) = event_signal(E, e_specs(:,e_s), dur);
            end
            
            %% Estimate FIR timecourses (marsbar manual) and fitted responses (CP)
            fprintf('Estimate event-related timecourses.\n')
            %% If your events have the same name across sessions, and you want to average across the events with the same name (if not, refer to manual for alternative code):
            % Get compound event types structure
            ets = event_types_named(E);
            n_event_types = length(ets);
            % Bin size in seconds for FIR
            bin_size = tr(E);
            % Number of FIR time bins to cover length of FIR
            bin_no = fir_length / bin_size;
            % Options - here �single� FIR model, return estimated % signal change
            opts = struct('single', 1, 'percent', 1);
            for e_t = 1:n_event_types
                fir_tc(:, e_t) = event_fitted_fir(E, ets(e_t).e_spec, bin_size, bin_no, opts);
                fitted_tc(:,e_t) = event_fitted(E,ets(e_t).e_spec,e_length);
            end
            
            %% Write tables (CP)
            fprintf('Write results to tables.\n')
            if pos == 0
                pos = pos+1;
                
                % contrast file
                % prepare header for contrast output file
                [nr_cons,~]=size(marsS.con);
                fprintf(con_file,'%s','subject');
                for j=1:nr_cons
                    fprintf(con_file,'\t%s',marsS.rows{1,j}.name); % print labels in header
                end
                % prepare variable coding output format of data
                con_format=['\n%s',repmat('\t%6.4f',1,nr_cons)];
                
                % percent signal change file
                % prepare header for contrast output file
                [~,nr_betas]=size(e_names);
                fprintf(psc_file,'%s','subject');
                for j=1:nr_betas
                    fprintf(psc_file,'\t%s',e_names{j}); % print labels in header
                end
                % prepare variable coding output format of data
                psc_format=['\n%s',repmat('\t%6.4f',1,nr_betas)];
                
            end
            
            % print contrast values to output file
            fprintf(con_file,con_format,subdirs(sub).name(1:dum1-1),marsS.con);
            
            % print psc values to output file
            fprintf(psc_file,psc_format,subdirs(sub).name(1:dum1-1),pct_ev);
            
            % write timecourse data to dummy matrix
            for e_t=1:n_event_types
                fir_tc_dum{pos}.tcs(e_t,:)=fir_tc(:,e_t); % timecourse data
                fir_tc_dum{pos}.conditions{e_t}=ets(e_t).name; % get condition label
                fir_tc_dum{pos}.vp=subdirs(sub).name;
                fitted_tc_dum{pos}.tcs(e_t,:)=fitted_tc(:,e_t); % timecourse data
                fitted_tc_dum{pos}.conditions{e_t}=ets(e_t).name; % get condition label
                fitted_tc_dum{pos}.vp=subdirs(sub).name;
            end
            
            clear('D','R','Y','E','xCon','b','marsS','bin_size','opts','fir_tc','fitted_tc','pct_ev','nr_cons','nr_betas');

            % prepare header for timecourse output files
            fprintf(fir_tc_file,'%s\t%s','subject','condition_label');
            for j=1:bin_no
                fprintf(fir_tc_file,'\t%s',int2str(j)); % print labels to header
            end
            
            fprintf(fitted_tc_file,'%s\t%s','subject','condition_label');
            [~,n_lines] = size(fitted_tc_dum{1}.tcs); % look up nr of bins of first subject
            for j=1:n_lines
                fprintf(fitted_tc_file,'\t%s',int2str(j));
            end
            
            % prepare variables coding output format of data
            fir_tc_format=['\n%s\t%s',repmat('\t%6.4f',1,bin_no)];
            fitted_tc_format=['\n%s\t%s',repmat('\t%6.4f',1,n_lines)];
            
            % print timecourse data to output file
            for e_t=1:n_event_types
                for j=1:pos
                    fprintf(fir_tc_file,fir_tc_format,subdirs(sub).name(1:dum1-1),fir_tc_dum{j}.conditions{e_t},fir_tc_dum{j}.tcs(e_t,:));
                    fprintf(fitted_tc_file,fitted_tc_format,subdirs(sub).name(1:dum1-1),fitted_tc_dum{j}.conditions{e_t},fitted_tc_dum{j}.tcs(e_t,:));
                end
            end
            
        end
        
        %         save([roi_list(roi,:) '_' suffix '_fir_timecourse'],'fir_tc_dum');
        %         save([roi_list(roi,:) '_' suffix '_fitted_timecourse'],'fitted_tc_dum');
        
        clear('ets','tc_dum','n_event_types');
    end
    fclose(con_file);
    fclose(psc_file);
    fclose(fir_tc_file);
    fclose(fitted_tc_file);
end
rmpath(genpath(folder));