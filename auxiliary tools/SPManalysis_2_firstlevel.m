%% Estimate spm firstlevel with spm12
% Christian Paret, ZI-Mannheim, 2022

% Requires dataset according to BIDS standard and derivatives output from fmriPrep (tested with version 20.2.0)
% Requires dataset-table according to RaBIDS standard
%
% You should provide a (unzipped) smoothed, preprocessed image in the derivatives directory. Use the "SPManalysis_1_smooth_images.m" program (see RaBIDS\auxiliary tools); it will automatically provide the images in the needed naming format.
%
% Download program to dataset/code directory to run
%%

% Change Log
% 2021/07/29: to identify tsv-file including regressors, search word "confounds" was changed to wild card for  downwards compatibility with fMRIPrep versions before 20.2.0
% 2021/08/09: definition of variable minImagescol added (bug remove)
% 2022/01/21: spm-firstlevel job: implemented check whether average of beta-weights for contrast specification is different from 0. In the future this code could be changed to return more specific information about consistency of SOTS and conditions_TaskID.xlsx-file
% 2022/04/12: ObjectType ContrastType added: can be eoi or tcon

clear
clc

fprintf('------------------------------------------------------------\n------------- Running firstlevel SPM analysis. -------------\n------------------------------------------------------------\n\nNOTE: Make sure that derivative images were unzipped and smoothed before running this program.\n      To unzip and smooth consider using programs unzip_niftis.m and SPManalysis_1_smooth_images.m that are available in the RaBIDS auxiliary tools directory.\n\n');
%% User input required: Provide information about your imaging data here:
smkernel = '6'; % smoothing kernel of functional images
TR = 2; % Repetition Time
slices = 36; % number of slices per volume
use_trimmed = true; % use trimmed session data, if it exists? If set true, the program will trim session data according to definition with program find_spikes_and_trim.m

%% User input required: Select options for nuisance regression:
get_realign = true; % 6 realignment regressors (3 translation, 3 rotation parameters)
get_outlier = false; % motion outlier "spike" regressors

%% Read data from datasheet
data = readtable('datasheet.xlsx','ReadRowNames',true,'PreserveVariableNames',true,'NumHeaderLines',0);
userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));
DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};
datasetd = [data_analysis_path,filesep,'dataset'];
first_imageline = find(strcmp(data.Properties.RowNames,'first image'));
minImagescol = find(strcmp(data.Properties.VariableNames,'MinImages'));

if contains(data{first_imageline,userInputcol}{:},'y')
    first_image = 1; % X=1:[first_image-1] images are skipped in create sots step and were deleted from the dicom-directory (the delete-mode will be removed in a future RaBIDS version and should not be used
    fimgs_exist = 0;
else
    first_image = data{first_imageline,minImagescol}; % X=1:[first_image-1] images are skipped in create sots step
    fimgs_exist = 1;
end

%% Identify derivative directory to work on
derivd = dir([data_analysis_path,filesep,'derivatives',filesep,'fmriprep*']);

fprintf(['Found ',num2str(length(derivd)),' derivative directories:\n'])

for i=1:length(derivd)
    fprintf([num2str(i),') ',derivd(i).name,'\n'])
end

try
    [c,~] = size(derivd);
catch
    fprintf('Cannot read derivatives directory.\n')
    return
end

if c > 1
    selectd = input('Enter number of derivative directory (e.g. 1) to work on and press enter.\n');
else
    selectd = 1;
end
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

%% Task to work on
reqtask = input('What task to work on? Enter a single task name (e.g. ''faces'').\n');

%% Define contrasts
condfile = ['conditions_',reqtask,'.xlsx'];
condata = readtable(fullfile(datasetd,'code',condfile),'ReadRowNames',true,'ReadVariableNames',true,'NumHeaderLines',0);
conditionlines = find(contains(condata.Properties.RowNames,'Condition'));

Namecol = find(strcmp(condata.Properties.VariableNames,'Name'));
ContrastTypecol = find(strcmp(condata.Properties.VariableNames,'ContrastType'));
ContrastPlus1col = find(strcmp(condata.Properties.VariableNames,'ContrastPlus1'));
ContrastMinus1col = find(strcmp(condata.Properties.VariableNames,'ContrastMinus1'));

contrastlines = find(contains(condata.Properties.RowNames,'Contrast'));

%% Go through subjects and calculate firstlevel
allsubs = dir([derivd,filesep,'sub-*']);

for subject = 1:length(allsubs)
    if allsubs(subject).isdir
        allses = dir(fullfile(derivd,allsubs(subject).name,'ses-*'));
        
        for session = 1:length(allses)
            if allses(session).isdir
                fprintf(['\nProcessing ',allsubs(subject).name,' ',allses(session).name,'.\n']);
                
                %% Define paths and do couple of checks before creating firstlevel model
                
                % Is there a Stimulus Onset Times file defining the stimulus protocol?
                multicondf = [allsubs(subject).name,'_',allses(session).name,'_task-',reqtask,'_multicond.mat'];
                multicondp = fullfile(datasetd,allsubs(subject).name,allses(session).name,'func');
                if isfile(fullfile(multicondp,multicondf))
                    multicond_ok = 1;
                else
                    fprintf('    SOTs file not found. --> SKIP SESSION.\n');
                    multicond_ok = 0;
                end
                
                % Does a smoothed derivative image exist?
                derivname = [allsubs(subject).name,'_',allses(session).name,'_task-',reqtask,'_*desc-preproc_desc-s',smkernel,'_bold.nii'];
                derivp = fullfile(derivd,allsubs(subject).name,allses(session).name,'func');
                derivnii = dir(fullfile(derivp,derivname));
                fprintf(['    Expecting derivative image ',derivnii.name,'.\n']);
                if length(derivnii) == 1
                    fprintf('        Found derivative image.\n');
                    deriv_ok = 1;
                elseif isempty(derivnii)
                    fprintf('        Derivative image not found. --> SKIP SESSION.\n');
                    deriv_ok = 0;
                else
                    fprintf('        Too many matching derivative images. --> SKIP SESSION.\n');
                    deriv_ok = 0;
                end
                
                % Load nuisance regressors
                % Routine below looks through all existing confounds timeseries files to pick out the right one
                % A txt-file containing the regressors will be saved to the firstlevel directory
                conf_timeseriesf = dir(fullfile(derivp,'sub-*-confounds_*.tsv'));
                if ~isempty(conf_timeseriesf)
                    for funcs = 1:length(conf_timeseriesf)
                        pos1 = strfind(conf_timeseriesf(funcs).name,'task');
                        pos2 = strfind(conf_timeseriesf(funcs).name,'_desc-confounds');
                        taskid = conf_timeseriesf(funcs).name(pos1+5:pos2-1);
                        
                        if strcmp(reqtask,taskid)
                            conf_table = readtable(fullfile(derivp,conf_timeseriesf(funcs).name),'FileType','text','ReadRowNames',false,'PreserveVariableNames',true,'NumHeaderLines',0);
                            
                            R = []; % Fill R with nuisance regressors
                            if get_realign
                                % realignment regressors
                                R = [conf_table.trans_x, conf_table.trans_y, conf_table.trans_z, conf_table.rot_x, conf_table.rot_y, conf_table.rot_z];
                            end
                            
                            if get_outlier
                                % motion outlier "spike" regressors
                                outliercol = find(contains(conf_table.Properties.VariableNames,'motion_outlier'));
                                for col = 1:length(outliercol) % this loop iterates through the motion_outlier colomns and concatenates R
                                    eval(['R = [R, conf_table.',conf_table.Properties.VariableNames{outliercol(col)},'];']);
                                end
                            end
                            %                         save(fullfile(sesd,'func',[allsubs(sub).name,'_',allses(ses).name,'_task-',taskid,'_desc-confounds',suffix,'_timeseries.mat']),'R')
                            fprintf(['    Confounds found for ',conf_timeseriesf(funcs).name,'.\n'])
                        end
                        multireg_ok = true;
                    end
                else
                    fprintf('     No confound file found. --> SKIP SESSION.\n');
                    multireg_ok = false;
                end
                
                % Check for trimsession file and define directory to write first-level SPM.mat
                firstleveld = [allsubs(subject).name,'_',allses(session).name,'_task-',reqtask];
                firstlevelp = fullfile(data_analysis_path,'spm_analysis','firstlevel',allses(session).name,['task-',reqtask],'taskrelated_activity'); % recommend suffix "taskrelated_standard" for future analyses
                trimf = [allsubs(subject).name,'_',allses(session).name,'_task-',reqtask,'_trimsession.mat'];
                trim = false; % set flag false as standard. Do now couple of checks to give detailed user feedback:
                
                if exist([derivp,filesep,trimf],'file') && use_trimmed
                    fprintf('    Found trimsession file.\n    Proceed with trimmed session data.\n')
                    trim = true;
                    firstlevelp = fullfile(firstlevelp,'trimmed_original',firstleveld);
                elseif exist([derivp,filesep,trimf],'file') && ~use_trimmed
                    fprintf('    Trimsession file exists, but use_trimmed is false.\n    Proceed with full session data.\n')
                    firstlevelp = fullfile(firstlevelp,'original',firstleveld);
                elseif ~exist([derivp,filesep,trimf],'file') && use_trimmed
                    fprintf('    Use_trimmed is true but no trimsession file was found.\n    Proceed with full session data.\n')
                    firstlevelp = fullfile(firstlevelp,'original',firstleveld);
                elseif ~exist([derivp,filesep,trimf],'file') && ~use_trimmed
                    fprintf('    Proceed with full session data.\n')
                    firstlevelp = fullfile(firstlevelp,'original',firstleveld);
                end
                
                % Check for firstlevel directory and if it does not exist: create it
                if exist(firstlevelp,'dir')
                    % Is there already an SPM.mat?
                    if ~isfile(fullfile(firstlevelp,'SPM.mat')) 
                        spm_ok = 1;
                    else
                        fprintf('\n --> Model appears to be estimated already, as an SPM.mat exists in the target directory. --> SKIP SESSION.\n');
                        spm_ok = 0;
                    end
                else
                    mkdir(firstlevelp)
                    spm_ok = 1;
                end
                
                if spm_ok && multireg_ok && multicond_ok && deriv_ok
                    fprintf('    \nAll model checks done. Now checking contrasts, if there are any defined...\n')
                    
                    %% Make and execute batch
                    % Define spm12 model
                    matlabbatch{1}.spm.stats.fmri_spec.dir = {firstlevelp};
                    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % do not change unless no conflict is expected with trimmed session data
                    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = slices;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = round(slices/2); % Input reference slice to which images were aligned to during slice-timing step. fMRIPrep aligns to middle slice if corresponding metadata information is available.
                    
                    getscans = cellstr(spm_select('Expand',fullfile(derivd,allsubs(subject).name,allses(session).name,'func',derivnii.name)));
                    
                    % if scan was trimmed: load trimsession file and adjust Stimulus Onset Times (SOTs) if necessary
                    % Note: Estimation of spm may crash if trimming leads us to
                    % discard all stimuli of a condition - this has not
                    % been tested and my require to clear empty cells after
                    % adjustment of stimulus timing below
                    if trim
                        load(fullfile(derivp,trimf)) % 1-by-2 vector with information of first/last scan to use; see script "find_spikes_and_trim_session.m"
                        firstscan = trimses(1);
                        lastscan = trimses(2);
                        
                        trimdelay = TR * (firstscan - first_image); % onset timing needs to be corrected by trimdelay; first_image denotes the reference scan (i.e., t0) for onsets created with RaBIDS_Create_SOTs
                        
                        SOTS = load(fullfile(multicondp,multicondf)); % Load stimulus onset times
                        names = SOTS.names;
                        
                        % if scans are trimmed in the beginning we need to correct onset of stimulus timing
                        if firstscan > first_image % i.e., trimdelay is positive; trimmed images are behind the images that have to be discarded
                            trim1_corrected_onsets = cell(1,length(SOTS.onsets)); % initiate a new onset matrix to fill
                            trim1_corrected_durations = cell(1,length(SOTS.durations)); % initiate a new onset matrix to fill
                            for k = 1:length(SOTS.onsets)
                                onsets{k} = SOTS.onsets{k} - trimdelay; % adjust onsets to cut "nose" of scans
                                for l = 1:length(onsets{k})
                                    if onsets{k}(l) < 0 % trimming removes sampled brain activation that relates to onset time of this event
                                        if onsets{k}(l) + SOTS.durations{k}(l) < 0 % the stimulus ended before the first volume that we include in our analysis
                                            fprintf(['        CAVE: Stimulus #',num2str(l),' of condition ',SOTS.names{k},' was removed by trimming.\n'])
                                        else
                                            trim1_corrected_onsets{k} = [trim1_corrected_onsets{k} 0]; % adjust onset time...
                                            trim1_corrected_durations{k} = [trim1_corrected_durations{k} onsets{k}(l)+SOTS.durations{k}(l)]; % ... and duration
                                            fprintf(['        CAVE: Onset time of stimulus #',num2str(l),' of condition ',SOTS.names{k},' was adjusted because of trimming.\n'])
                                        end
                                    else % change nothing and write to matrix
                                        trim1_corrected_onsets{k} = [trim1_corrected_onsets{k} onsets{k}(l)];
                                        trim1_corrected_durations{k} = [trim1_corrected_durations{k} SOTS.durations{k}(l)];
                                    end
                                end
                            end
                            
                            % save adjusted onsets and durations matrices to workspace
                            onsets = trim1_corrected_onsets;
                            durations = trim1_corrected_durations;
                            
                        elseif firstscan <= first_image % if trimmed scans fall under the volumes to be discarded (to account for T1-effects in beginning of session) ...
                            firstscan = first_image; % ... use the image identified by first_image as first scan in the model; first_image denotes the reference scan (i.e., t0) for onsets created with RaBIDS_Create_SOTs
                            
                            % save original onsets and durations matrices to workspace
                            onsets = SOTS.onsets;
                            durations = SOTS.durations;
                        end
                        
                        if lastscan < length(getscans) % check whether volumes are trimmed in the end of the scan
                            fprintf('          CAVE: Adjustment of SOTs for tail-trimmed scans is not sufficiently tested.\n')
                            truelength = (lastscan * TR) - trimdelay; % this is the length of the scan after correction for trimming in the beginning
                            trim2_corrected_onsets = cell(1,length(SOTS.onsets)); % initiate a new onset matrix to fill
                            trim2_corrected_durations = cell(1,length(SOTS.durations)); % initiate a new onset matrix to fill
                            for k = 1:length(onsets)
                                for l = 1:length(onsets{k})
                                    if onsets{k}(l) + SOTS.durations{k}(l) > truelength % duration of stimulus event extends available volumes
                                        if onsets{k}(l) > truelength % onset of stimulus is presented during scans that are to be discarded through trimming
                                            fprintf(['        CAVE: Stimulus #',num2str(l),' of condition ',SOTS.names{k},' was removed by trimming.\n'])
                                        else
                                            trim2_corrected_onsets{k} = [trim2_corrected_onsets{k} onsets{k}(l)]; % keep onset as is ...
                                            trim2_corrected_durations{k} = [trim2_corrected_onsets{k} SOTS.durations{k}(l)-trimdelay]; % ... and adjust duration
                                            fprintf(['        CAVE: Duration of stimulus #',num2str(l),' of condition ',SOTS.names{k},' was adjusted because of trimming.\n'])
                                        end
                                    else % change nothing and write to matrix
                                        trim2_corrected_onsets{k} = [trim2_corrected_onsets{k} onsets{k}(l)];
                                        trim2_corrected_durations{k} = [trim2_corrected_durations{k} SOTS.durations{k}(l)];
                                    end
                                end
                            end
                            
                            % save adjusted onsets and durations matrices to workspace
                            onsets = trim2_corrected_onsets;
                            durations = trim2_corrected_durations;                            
                        end
                        
                        multicondf = [allsubs(subject).name,'_',allses(session).name,'_task-',reqtask,'_trimmed_multicond.mat'];
                        save(fullfile(multicondp,multicondf),'names','onsets','durations') % save a multicond-file with trimmed SOTs
                        
                    else
                        firstscan = first_image;
                        lastscan = length(getscans);
                    end
                    
                    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = getscans(firstscan:lastscan);                     
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {fullfile(multicondp,multicondf)};
                    
                    % write nuisance text file to firstlevel directory; this is a workaround
                    writematrix(R(firstscan:lastscan,:),fullfile(firstlevelp,'nuisance_regressors.txt'),'Delimiter','tab');
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(firstlevelp,'nuisance_regressors.txt')};
                    
                    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
                    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
                    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
                    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
                    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.3;
                    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
                    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
                    
                    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
                    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
                    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
                    
                    multiconddata = load(matlabbatch{1}.spm.stats.fmri_spec.sess.multi{:});
                    
                    % get weights and contrast type information, check validity of contrasts                    
                    for i = 1:length(contrastlines)
                        contype = condata{contrastlines(i),ContrastTypecol}{:};
                            
                        if strcmp(contype,'eoi')
                            fprintf(['    \n    Contrast ',condata{contrastlines(i),Namecol}{:},' is eoi F-contrast. Check below whether inclusion/exclusion of conditions is as expected.\n']);
                            contype = 'fcon';
                            
                            % get conditions to include in eoi
                            dum1 = condata{contrastlines(i),ContrastPlus1col}{:};
                            condstr = dum1(~isspace(dum1));
                            separator = [0 strfind(condstr,';')]; % conditions are separated by semicolon
                            if length(separator) > 1
                                for j = 1:length(separator)-1
                                    eoicond{j} = condstr(separator(j)+1:separator(j+1)-1);
                                end
                                eoicond{j+1} = condstr(separator(j+1)+1:end);
                            else
                                eoicond{1} = condstr;
                            end
                            
                            % make weight-matrix with as many rows as there are conditions, and as many colomns as there are eois
                            eoimatrix = eye(length(conditionlines));
                            
                            includecon = [];
                            for j = 1:length(conditionlines)
                                checkcond = condata{conditionlines(j),Namecol}{:};
                                fprintf(['      Found condition ''',checkcond,'''.\n'])
                                for k = 1:length(eoicond)
                                    if strcmp(checkcond,eoicond{k})
                                        includecon = [includecon j];
                                        fprintf(['        --> include ''',checkcond,''' in eoi f-contrast.\n'])
                                    end
                                end
                            end
                            
                            if includecon
                                eoimatrix = eoimatrix(includecon,:);
                            else
                                fprintf('      --> No conditions to include in f-contrast. Check whether conditions are entered correctly.\nError #17.\n')
                                return
                            end
                                                        
                            weights = eoimatrix;
                            
                        elseif strcmp(contype,'tcon')
                            fprintf(['    \n    Contrast ',condata{contrastlines(i),Namecol}{:},' is a t-contrast.\n']); 
                        
                            weights = [];
                            for j = 1:length(multiconddata.names)
                                if strcmp(condata{contrastlines(i),ContrastPlus1col}{:},multiconddata.names{j})
                                    weights = [weights 1];
                                elseif strcmp(condata{contrastlines(i),ContrastMinus1col}{:},multiconddata.names{j})
                                    weights = [weights -1];
                                elseif sum(weights)~=0
                                    fprintf(['    -------- CAVE -------\n        Weights of contrast ',condata{contrastlines(i),Namecol}{:},' do not average to 0.\n        If this is unexpected it is recommended to check ',condfile,' for consistency with stimulus onset times (sots)!    \n------------------------\n']);
                                else
                                    weights = [weights 0];
                                end
                            end

                            if ~any(weights)
                                fprintf(['    Contrast ',condata{contrastlines(i),Namecol}{:},' is invalid.\n']);
                                return
                            end
                            
                        else
                            fprintf(['    ContrastType of contrast ',condata{contrastlines(i),Namecol}{:},' is invalid. Can be eoi or tcon.\n']);
                            return
                        end
                        
                        conmat = mat2str(weights);
                        
                        % fill batch with contrast information
                        eval(['matlabbatch{3}.spm.stats.con.consess{',num2str(i),'}.',contype,'.name = ''',condata{contrastlines(i),Namecol}{:},''';']);
                        eval(['matlabbatch{3}.spm.stats.con.consess{',num2str(i),'}.',contype,'.weights = ',conmat,';']);
                        eval(['matlabbatch{3}.spm.stats.con.consess{',num2str(i),'}.',contype,'.sessrep = ''bothsc'';']); % can be 'none' in case of single-session data
                    end
                    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
                    matlabbatch{3}.spm.stats.con.delete = 1;
                    
                    %% Estimate model
                    fprintf('\nEstimate SPM model...\n');
                    spm_jobman('run', matlabbatch);
                    fprintf('Model estimation successful!\nResults can be found here: %s.\n',firstlevelp);
                    clear matlabbatch
                else
                    fprintf('\n    At least one of the checks above found an issue. Skipping session.\n')
                end
            end
        end
    end
end

fprintf('\nEnd.\n')