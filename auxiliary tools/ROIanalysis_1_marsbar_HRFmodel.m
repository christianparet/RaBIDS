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
bin_size = 2; % define bin_size in seconds. 
bufferseconds = 20; % add time (s) to sample BOLD signal following stimulus presentation
writeheadercon = 1;
writeheaderpsc = 1;
firstleveldir = uigetdir(pwd,'Select firstlevel directory including the subject directories'); % directory containing child directories with each child containing an estimated firstlevel SPM.mat
subdirs = dir([firstleveldir,filesep,'sub*']);
maskdir = uigetdir(pwd,'Select directory containing binary brain mask files (*roi.mat)');

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
    filename = fullfile(firstleveldir,'readme.txt');
    readme_file = fopen(filename,'w');
    
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
            
            spm_name = fullfile(firstleveldir,subdirs(sub).name,'SPM.mat');
            
            if exist(spm_name,'file')
                
                spm_data = load(spm_name); % load spm to get stimulus timings
                
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
                fprintf('Extract percent signal change and estimate event-related timecourses.\n')
                
                try
                    [e_specs, e_names] = event_specs(D);
                    n_events = size(e_specs, 2);
                    
%                     bin_size = tr(E)*2; % Bin size in seconds for FIR
                    opts = struct('single', 1, 'percent', 1); % Options - here �single� FIR model, return estimated % signal change
                    
                    for d_s = 1:n_events
                        alldur(d_s) = mean(spm_data.SPM.Sess.U(:,1).dur); % Caution: assuming that all events of a condition have similar duration
                    end
                    
                    dur = max(alldur)+bufferseconds; % set duration to sample timecourses to maximum length needed. 
                    
                    notcdata = []; % initiate matrix to collect names of events (i.e. tasks) for which timecourses cannot be estimated
                    
                    for e_s = 1:n_events
                        
                        pct_ev(e_s) = event_signal(E, e_specs(:,e_s), dur); % Extract percent signal change
                        
                        try
                            bin_no = round(dur / bin_size); % Number of FIR time bins to cover length of FIR
                            fir_tc(:, e_s) = event_fitted_fir(E, e_specs(:,e_s), bin_size, bin_no, opts); % FIR timecourse
                            fitted_tc(:,e_s) = event_fitted(E,e_specs(:,e_s),dur); % Fitted timecourse
                        catch
                            fprintf(['   An error occured: cannot get FIR and fitted timecourse for ',e_names{e_s},'. This may happen when a task event is interrupted e.g. because of disrupted scan.\n']);
                            notcdata = [notcdata; e_names(e_s)];
                        end
                    end
                    psc_success = true;
                catch
                    fprintf('   Failed. This can happen because no events have been modeled.\n')
                    psc_success = false;
                end
                
                %% Write tables (CP)
                fprintf('Write results to tables.\n')
                
                if writeheadercon == 1
                    writeheadercon = 0;
                    % contrast file
                    % prepare header for contrast output file
                    [nr_cons,~]=size(marsS.con);
                    fprintf(con_file,'%s','subject');
                    for j=1:nr_cons
                        fprintf(con_file,'\t%s',marsS.rows{1,j}.name); % print labels in header
                    end
                    % prepare variable coding output format of data
                    con_format=['\n%s',repmat('\t%6.4f',1,nr_cons)];
                end
                
                % print contrast values to output file
                fprintf(con_file,con_format,subdirs(sub).name(1:dum1-1),marsS.con);
                
                if psc_success                                       
                    % write timecourse data to dummy matrix
                    for e_t=1:n_events
                        if ~any(strcmp(notcdata,e_names(e_t)))
                            fir_tc_dum.tcs(e_t,:)=fir_tc(:,e_t); % timecourse data
                            fir_tc_dum.conditions(e_t)=e_names(e_t); % get condition label
                            %                             fir_tc_dum.vp=subdirs(sub).name(1:dum1-1);
                            fitted_tc_dum.tcs(e_t,:)=fitted_tc(:,e_t); % timecourse data
                            fitted_tc_dum.conditions(e_t)=e_names(e_t); % get condition label
                            %                             fitted_tc_dum.vp=subdirs(sub).name(1:dum1-1);
                        end
                    end
                end
                
                clear('D','R','Y','E','xCon','b','marsS','opts','fir_tc','fitted_tc','nr_cons','nr_betas');
                
                if psc_success
                    
                    if writeheaderpsc == 1
                        writeheaderpsc=0;
                        
                        % percent signal change file
                        % prepare header for contrast output file
                        [~,nr_betas]=size(e_names);
                        fprintf(psc_file,'%s','subject');
                        for j=1:nr_betas
                            fprintf(psc_file,'\t%s',e_names{j}); % print labels in header
                        end
                        % prepare variable coding output format of data
                        psc_format=['\n%s',repmat('\t%6.4f',1,nr_betas)];
                        
                        % prepare header for timecourse output files
                        fprintf(fir_tc_file,'%s\t%s','subject','condition_label');
                        for j=1:bin_no
                            fprintf(fir_tc_file,'\t%s',int2str(j)); % print labels to header
                        end
                        
                        fprintf(fitted_tc_file,'%s\t%s','subject','condition_label');
                        [~,n_lines] = size(fitted_tc_dum.tcs); % look up nr of bins of first subject
                        for j=1:n_lines
                            fprintf(fitted_tc_file,'\t%s',int2str(j));
                        end
                        
                        % prepare variables coding output format of data
                        fir_tc_format=['\n%s\t%s',repmat('\t%6.4f',1,bin_no)];
                        fitted_tc_format=['\n%s\t%s',repmat('\t%6.4f',1,n_lines)];
                    end
                    
                    % print psc values to output file
                    fprintf(psc_file,psc_format,subdirs(sub).name(1:dum1-1),pct_ev);
                    
                    % print timecourse data to output file
                    for e_t=1:n_events
                        if ~any(strcmp(notcdata,e_names(e_t)))
                            fprintf(fir_tc_file,fir_tc_format,subdirs(sub).name(1:dum1-1),fir_tc_dum.conditions{e_t},fir_tc_dum.tcs(e_t,:));
                            fprintf(fitted_tc_file,fitted_tc_format,subdirs(sub).name(1:dum1-1),fitted_tc_dum.conditions{e_t},fitted_tc_dum.tcs(e_t,:));
                        end
                    end
                    
                    clear('fir_tc_dum', 'fitted_tc_dum','pct_ev')
                else
                    fprintf('   No BOLD percent signal change and no time course data available.\n')
                    
                end
            else
                fprintf('   No SPM.mat found for this subject.\n')
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
    
    title = 'Data in ROI txt files produced with marsbar-0.4.4';
    ref_lit = 'Matthew Brett, Jean-Luc Anton, Romain Valabregue, Jean-Baptiste Poline. Region of interest analysis using an SPM toolbox [abstract] Presented at the 8th International Conference on Functional Mapping of the Human Brain, June 2-6, 2002, Sendai, Japan. Available on CD-ROM in NeuroImage, Vol 16, No 2.';
    if psc_success
        ref_web = 'ROI-percent-BOLD-signal-change txt file: To learn how percent signal change was calculated see http://marsbar.sourceforge.net/index.html';
        fir_met = ['ROI-FIR-timecourse txt file: Timecourse was modeled with a Finite Impulse Response (FIR) function, start is stimulus onset, with bin size ',int2str(bin_size),' s (per unit on x-axis) and duration of ',int2str(dur),' s.'];
        fit_met = ['ROI-fitted-timecourse txt file: Timecourse was modeled with the HRF, start is stimulus onset, with duration of ',int2str(dur),' s.'];
    end
    ref_rabids = 'Extraction of data was run with an in-house matlab script called ROIanalysis_1_extract.m, available via https://github.com/christianparet/RaBIDS in the auxiliary tools directory.';
    
    if psc_success
        fprintf(readme_file,'%s\n\n%s\n%s\n%s\n\n%s\n%s\n%s',title, ref_web,fir_met,fit_met,'References:',ref_lit,ref_rabids);
    else
        fprintf(readme_file,'%s\n\n%s\n%s\n%s',title, 'References:',ref_lit,ref_rabids);
    end
    fclose(readme_file);
end
rmpath(genpath(folder));