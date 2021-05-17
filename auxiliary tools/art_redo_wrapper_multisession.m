
%% art_redo_wrapper
% art_redo program, modified by Christian Paret, ZI-Mannheim, 2021/04/07

% This program runs art_redo from the ArtRepair v5b toolbox over data; it
% re-estimates an existing firstlevel SPM model using repaired data
% For each original and repaired SPM model it calls art_summary to estiamte
% Global Quality metrics and plot Residuals.

% To use it you first need to repair data; use e.g. the Repair_volumes.m program, which is a wrapper for art_global and is available in RaBIDS/auxiliary tools on my github

clear

%%
% Specify name of directories to look for data and to save it. The data structure arises from the RaBIDS pipeline and the Repair_volumes.m program
spm_analysis_dir = 'Y:\Projects\EFPTest\Data_analysis\spm_analysis';
threednii =  'repaired_preprocessed_3dnii'; % path to repaired data
sesdir = 'ses-post'; % the session you want to process
taskdir = 'task-AllAvailable'; % the name of the firstlevel task-directory (can be different than single-task names/taskID, e.g. in multi-session design)

taskspecdir = 'efpmodel'; % defines subdirectory of firstlevel task-directory, where the original SPM can be found

conimage = 'con_0004.nii'; % define the con image to compute quality indices
backupcon = 'con_0001.nii'; % this option is introduced as workaround for subjects that do not have con_0004.nii (because of missing Session data)

% subdirectories of threednii directory are structured:
% <spm_analysis_dir>/<threednii>/<sesdir>/<subjectID>/<taskdir>/...

trim = true; % if repaired + trimmed 3d-nii data exist set trim=true. trim=false will process complete session data only.
%%

clc

% function art_redo
% FUNCTION art_redo   v3.2
%
% Applies repairs and deweighting to an existing SPM design. This program
% asks to replace the images in each session with repaired images, and
% applies deweighting in estimation according to the art_deweighted.txt file.
% It assumes user has repaired images with art_global to create repaired
% images (with prefix v). Each session should be repaired separately.
%
% Deweighting is applied as a fixed value of 0.01 on all repaired images,
% essentially knocking them out of the estimate, but maintaining the size
% of the design matrix. Deweighted scans show up as dark bars on the
% SPM design matrix display. User can compare before and after results
% with Global Quality program.
%
% INPUT by GUI
%    Specify the SPM.mat file.
%    Specify a new folder for repaired Results.
%    Specify the replacement images for each session
% OUTPUT
%    Writes a new SPM.mat in the designated new Results folder,
%       with deweight factors and the replacement images.
%    Runs SPM estimation and contrasts for this repaired design,
%       and writes beta, con, etc. images in new Results folder.
%
% All original SPM.mat and results files are preserved.
% BUGS: Hardly affects the estimates since repairs to do most of work.
%       Prevents the whitening process in SPM estimation.
%
% v3.1  May09 pkm
% v3.2 supports SPM12 Dec2014 pkm


% v2.2  changed SPM.xVi.form from 'i.i.d' to 'none'.   Jul 2007
%       scales to % using peak of regressor and contrast sum.
% v3.   Removed functionality to repair and run global quality.
%       assumes each file realigned and repaired separately. Mar09 pkm
% v3.1  fixed bugs (RMitchell): Deweight=1, spm_ver
% Paul Mazaika, Jan 2007
% Bug fix for output. A. Gonzalo

clear SPM scans;
% Configure while preserving old SPM versions
spmv = spm('Ver'); spm_ver = 'spm5';  % chooses spm_select to read vols
if (strcmp(spmv,'SPM2')) spm_ver = 'spm2'; end
if (strcmp(spmv,'SPM2') || strcmp(spmv,'SPM5')) spm_defaults;
else spm('Defaults','fmri'); end

addpath('C:\Program Files\spm12\toolbox\ArtRepair')
% List all folders in directory <curr dir>\regular and use art_redo on content, by CP 2021/04/06
allsubs = dir(fullfile(spm_analysis_dir,'firstlevel',sesdir,taskdir,taskspecdir,'regular','sub-*'));

subcount = 1;

Deweight = 1;

for sub = 1:length(allsubs)
    goahead = 1;
    fprintf(['\nProcessing ',allsubs(sub).name,'\n'])
    
    % find subject string and write it in varialbe subID
    pos = strfind(allsubs(sub).name,'_ses-');
    subID = allsubs(sub).name(1:pos-1);
    
    repfirstlevel = fullfile(spm_analysis_dir,'firstlevel',sesdir,taskdir,taskspecdir,'repaired'); % this is the directory where we write the repaired firstlevel model to
    origSPMpath = fullfile(spm_analysis_dir,'firstlevel',sesdir,taskdir,taskspecdir,'regular',allsubs(sub).name); % this is the directory containing the original model (SPM.mat) that art_redo should use to re-estimate with repaired volumes

    if exist(fullfile(origSPMpath,'SPM.mat'),'file')
        fprintf('Found original firstlevel SPM.\n')
        spm_available = 1;
        
        datafile_ID(subcount,:) = allsubs(sub).name;
        session_data(subcount,:) = 'complete';
        
        % Define scale factors to be passed to art_summary (I re-used code found in the art_summary function
        imgmask = fullfile(origSPMpath,'mask.nii');
        try
            ScaleFactors = art_percentscale(fullfile(origSPMpath,conimage),imgmask);
        catch
            ScaleFactors = art_percentscale(fullfile(origSPMpath,backupcon),imgmask);
            fprintf([conimage,' not found, use ',backupcon,' instead for global quality assessment.\n'])
        end
        peak_value = ScaleFactors(1);
        contrast_value = ScaleFactors(2);
        bmean = ScaleFactors(3);
        
        try
            [StdEstimationError_origSPM(subcount,:),MeanEstimationError_origSPM(subcount,:),ResMSerror_origSPM(subcount,:)] = art_summary(fullfile(origSPMpath,conimage),imgmask,origSPMpath,allsubs(sub).name,1,ScaleFactors);
        catch
            [StdEstimationError_origSPM(subcount,:),MeanEstimationError_origSPM(subcount,:),ResMSerror_origSPM(subcount,:)] = art_summary(fullfile(origSPMpath,backupcon),imgmask,origSPMpath,allsubs(sub).name,1,ScaleFactors);
        end
        h(1) = figure(1);
        h(2) = figure(2);
        savefig(h,fullfile(origSPMpath,'ResidualPlots.fig'))
        close(h)
        
    else
        fprintf('Original firstlevel SPM not found. Skip session.\n')
        goahead = 0;
        spm_available = 0;
    end
    
    if exist(fullfile(repfirstlevel,allsubs(sub).name),'dir')
        if exist(fullfile(repfirstlevel,allsubs(sub).name,conimage),'file')
            fprintf('Model has already been estimated with repaired data. Skip subject.\n')
            goahead = 0;
            
            % Define scale factors to be passed to art_summary (I re-used code found in the art_summary function
            imgmask = fullfile(origSPMpath,'mask.nii'); % use the same mask as was used for original SPM
            try
                ScaleFactors = art_percentscale(fullfile(repfirstlevel,allsubs(sub).name,conimage),imgmask);
            catch
                ScaleFactors = art_percentscale(fullfile(repfirstlevel,allsubs(sub).name,backupcon),imgmask);
                fprintf([conimage,' not found, use ',backupcon,' instead for global quality assessment.\n'])
            end
            peak_value = ScaleFactors(1);
            contrast_value = ScaleFactors(2);
            bmean = ScaleFactors(3);
            
            try
                [StdEstimationError_redoSPM(subcount,:),MeanEstimationError_redoSPM(subcount,:),ResMSerror_redoSPM(subcount,:)] = art_summary(fullfile(repfirstlevel,allsubs(sub).name,conimage),imgmask,fullfile(repfirstlevel,allsubs(sub).name),allsubs(sub).name,1,ScaleFactors);
            catch
                [StdEstimationError_redoSPM(subcount,:),MeanEstimationError_redoSPM(subcount,:),ResMSerror_redoSPM(subcount,:)] = art_summary(fullfile(repfirstlevel,allsubs(sub).name,backupcon),imgmask,fullfile(repfirstlevel,allsubs(sub).name),allsubs(sub).name,1,ScaleFactors);
            end
            h(1) = figure(1);
            h(2) = figure(2);
            savefig(h,fullfile(repfirstlevel,allsubs(sub).name,'ResidualPlots.fig'))
            close(h)
            
        end
        
    else
        mkdir(fullfile(repfirstlevel,allsubs(sub).name))
    end
    
    if spm_available % introduced by CP, 2021/04/06
        
        % Get the existing SPM file
        
        %     if strcmp(spm_ver,'spm5') % Commented out by CP, 2021/04/06
        %         origSPM   = spm_select(1,'mat','Select SPM design to re-estimate'); % Commented out by CP, 2021/04/06
        %         RepairResults = spm_select(1,'dir','Select folder for Repaired Results'); % Commented out by CP, 2021/04/06
        
        origSPM   = spm_select('FPList',origSPMpath,'SPM.mat'); % by CP, 2021/04/06
        RepairResults = spm_select('FPList',repfirstlevel,'dir',allsubs(sub).name); % by CP, 2021/04/06
        
        %     else % Commented out by CP, 2021/04/06
        %         origSPM = spm_get(1, 'SPM.mat', 'Select SPM design to re-estimate'); % Commented out by CP, 2021/04/06
        %         RepairResults = spm_get(-1,'*','Select folder for Repaired Results'); % Commented out by CP, 2021/04/06
        %     end % Commented out by CP, 2021/04/06
        
        dirSPM = fileparts(origSPM);
        %     cd(dirSPM); % Commented out by CP, 2021/04/06
        load(fullfile(origSPMpath,'SPM'))
        num_sess = size(SPM.Sess,2);  % size(SPM.nscan,2)
               
        [Nvol, ~] = size(SPM.xY.P); % write number of volumes of original SPM to Nvol
    end
    
    if goahead % Access code below to redo SPM model with repaired data. Otherwise skip it. %introduced by CP, 2021/04/06
        
        % SELECT A NEW FOLDER FOR REPAIRED RESULTS
        SPM.swd = RepairResults;
        
        % Get the repaired data
        AutoRepair = 0;
        %spm_input('!DeleteInputObj');
        if AutoRepair == 0
            % (Preferable if program automatically checked for necessary v-files.)
            
            %% Update each session with existing v* file data.
%             nwd=cd(fileparts(SPM.xY.P(1,:))); % fails if upper paths % commented out by CP, 2021/05/11
            
            lastimg = 0;
            for j=1:num_sess
                lastimg = lastimg + SPM.nscan(j); % lastimg gets the list number of the last image of this Session
                [nwd, derivn] = fileparts(SPM.xY.P(lastimg,:)); % write information on location of volumes to variables % by CP, 2021/05/11
                dum1 = strfind(derivn,'task-'); % identify the task that belongs to SPM's session j % by CP, 2021/05/11
                dum2 = strfind(derivn,'_space-'); % by CP, 2021/05/11
                taskID = derivn(dum1:dum2-1); % by CP, 2021/05/11
                
                sessdata = dir(fullfile(spm_analysis_dir,threednii,sesdir,subID,taskID,'*_session_data')); % by CP, 2021/05/11
                if length(sessdata) > 1 % by CP, 2021/05/11
                    fprintf('    CAVE: Repaired images task-directory has two sub-directories but must have only one!') % by CP, 2021/05/11
                    cd(sessdata) % by CP, 2021/05/11
                    return % by CP, 2021/05/11
                elseif length(sessdata) < 1 % by CP, 2021/05/11
                    fprintf('    Repaired data was not found.\n') % by CP, 2021/05/11
                    return % by CP, 2021/05/11
                end % by CP, 2021/05/11
                
                repdatadir = fullfile(spm_analysis_dir,threednii,sesdir,subID,taskID,sessdata.name); % directory containing the repaired 3d-nii data % by CP, 2021/05/11
                
                Pnew{j} = spm_select('FPList',repdatadir,'^v.*' ); % by CP, 2021/05/11
                
                %                 askstring = ['Select repaired images for session ' num2str(j) ]; % Commented out by CP, 2021/04/06
                %                 if strcmp(spm_ver,'spm5') % Commented out by CP, 2021/04/06
                %                     Pnew{j}   = spm_select(Inf,'image', askstring,[],dirSPM,'^v.*' ); % Commented out by CP, 2021/04/06
                %                 else % Commented out by CP, 2021/04/06
                %                     Pnew{j}   = spm_get(Inf,'v*img', askstring, nwd );  % for SPM2 % Commented out by CP, 2021/04/06
                %                 end % Commented out by CP, 2021/04/06
            end
            
            %% To concatenate character arrays defining volumes: account for variable array size (CP 2021/05/11)
            Narrays=[];
            for i = 1:num_sess
                Narrays = [Narrays,min(numel(Pnew{i}(1,:)))]; % first, we collect sizes of character arrays in Narrays
            end
                        
            scans = repmat(' ',Nvol,max(Narrays)); % inititate character matrix containing whitespaces of needed size
            
            % loop through tasks (SPM-language/BIDS language: Session/Task) and fill matrix scans with repaired volumes
            dum = 0;
            for i = 1:num_sess
                scans(dum+1:dum+SPM.nscan(i),1:length(Pnew{i}(1,:))) = Pnew{i};
                dum = dum+SPM.nscan(i);
            end
            clear Pnew
            
            if ~(size(SPM.xY.P,1) == size(scans,1)) 
                disp('ERROR: Must enter same number of new data points.')
                return; 
            else     % Assign the repaired data to the design matrix. 
                SPM.xY.P           = scans; 
            end
            %%
        end
        
        % Configure the design matrix with the new scans  (v2.2)
        if (strcmp(SPM.xVi.form, 'i.i.d'))  % changed name is used in spm_fmri_spm_ui.
            SPM.xVi.form = 'none';
        end
        SPM = spm_fmri_spm_ui(SPM);
        
        %=====================================================================
        %  LOGIC FOR DEWEIGHTING SCANS
        %-------------------------------------------------------------------------
        %    If SPM.xX.W is not provided, then SPM will default to starting
        % the analysis with SPM.xX.W = Identity Matrix, and will update it
        % for whitening. There is only one SPM.xX.W input that covers all
        % the sessions for a multi-session study.
        %    When SPM.xX.W is provided as input, then SPM uses the initial
        % weighting in the analysis and SPM does not do more whitening.
        %    When the ArtifactRepair program "repairs" the data, it writes the
        % files art_deweighted.txt and art_repaired.txt to the Images folder in
        % a single session study. For multiple sessions, we assume each session
        % is repaired separately. If those files exist, the logic below
        % initializes SPM.xX.W.
        %
        % Possible alternatives:
        %  -To only deweight the repaired files, read art_repaired.txt instead.
        %  -The omit scans logic in the 2005 version was here instead of deweight.
        
        if Deweight == 1
            % For multiple sessions, deweightlist is in EACH session folder.
            Vom = ones(1,size(scans,1));
            try
                firstimg = 1; % approach to identify number of scans per session introduced by CP, 2021/05/12; original ArtRepair code assumes equal number of scans per Session
                for isess = 1:num_sess
                    imagedir = fileparts(scans(firstimg,:)); % write information on location of volumes to variables % by CP, 2021/05/11
                    deweightlist = fullfile(imagedir,'art_deweighted.txt');
                    outindex = load(deweightlist);
                    Vom(outindex + firstimg -1) = 0.01;
                    outmsg = [num2str(length(outindex)) ' images will be deweighted in session ' num2str(isess)];
                    disp(outmsg)
                    firstimg = firstimg + SPM.nscan(isess); % lastimg gets the list number of the last image of this Session
                end
                SPM.xX.W = sparse(diag(Vom,0));
                disp('Applying Deweighting to GLM Estimatees')
            catch
                dummy = 1;  %  SPM.xX.W will follow SPM defaults.
                disp('WARNING: Could not locate a deweighting file for each session')
            end
        end
        
        % %================================
        % %  OLD LOGIC LOGIC FOR DEWEIGHTING SCANS
        % %--------------------------------
        % % This version assumed all sessions were realigned together.
        % deweightdir = fileparts(SPM.xY.P(1,:));
        % deweightlist = fullfile(deweightdir,'art_deweighted.txt');
        %
        % Deweight = 1;
        % if Deweight == 1
        %     try
        %         outindex = load(deweightlist);
        %         Vom = ones(1,size(SPM.xY.P,1));
        %         Vom(outindex) = 0.01;   % better than deweighting by 0.1
        %         SPM.xX.W = sparse(diag(Vom,0));
        %         disp('Applying Deweighting Factors to Images')
        %         outmsg = [num2str(length(outindex)) ' images will be deweighted'];
        %         disp(outmsg)
        %     catch
        %         dummy = 1;  %  SPM.xX.W will follow SPM defaults.
        %         disp('WARNING: Could not locate a deweighting file')
        %     end
        % end
        
        %-------------------------------------------------------------
        % Write the amended SPM.mat file with deweighting and new data included.
        disp('The SPM.mat file has been modified and saved to disk.');
        cd(SPM.swd);
        save SPM SPM;
        
        
        % Estimate parameters and contrasts
        
        % Save contrast definitions, before deletion by spm_spm.
        temp = SPM.xCon;
        consize = length(SPM.xCon);
        SPM = spm_spm(SPM);
        
        % Recover the contrasts and run them
        SPM.xCon = temp;  % explicitly sets up structure for SPM5, but fills in too much.
        for j = 1:consize
            spmtemp = spm_FcUtil('Set', temp(j).name, temp(j).STAT, 'c', temp(j).c, SPM.xX.xKXs);
            SPM.xCon(j) = spmtemp;
        end
        if consize > 0
            spm_contrasts(SPM);
        end
        
        disp('Done. New estimates have been created.');
        disp('Compare new and old estimates with Global Quality metric.');
        
        % Define scale factors to be passed to art_summary (I re-used code found in the art_summary function
        imgmask = fullfile(origSPMpath,'mask.nii'); % use the same mask as was used for original SPM
        
        try
            ScaleFactors = art_percentscale(fullfile(repfirstlevel,allsubs(sub).name,conimage),imgmask);
        catch
            ScaleFactors = art_percentscale(fullfile(repfirstlevel,allsubs(sub).name,backupcon),imgmask);
            fprintf([conimage,' not found, use ',backupcon,' instead for global quality assessment.\n'])
        end
        peak_value = ScaleFactors(1);
        contrast_value = ScaleFactors(2);
        bmean = ScaleFactors(3);
        
        try
            [StdEstimationError_redoSPM(subcount,:),MeanEstimationError_redoSPM(subcount,:),ResMSerror_redoSPM(subcount,:)] = art_summary(fullfile(repfirstlevel,allsubs(sub).name,conimage),imgmask,fullfile(repfirstlevel,allsubs(sub).name),allsubs(sub).name,1,ScaleFactors);
        catch
            [StdEstimationError_redoSPM(subcount,:),MeanEstimationError_redoSPM(subcount,:),ResMSerror_redoSPM(subcount,:)] = art_summary(fullfile(repfirstlevel,allsubs(sub).name,backupcon),imgmask,fullfile(repfirstlevel,allsubs(sub).name),allsubs(sub).name,1,ScaleFactors);
        end        
        h(1) = figure(1);
        h(2) = figure(2);
        savefig(h,fullfile(repfirstlevel,allsubs(sub).name,'ResidualPlots.fig'))
        close(h)
    end
    
    if spm_available % introduced by CP, 2021/06/04
        
        if Deweight == 1
            repSPM = load(fullfile(RepairResults,'SPM.mat'));
            
            % For multiple sessions, deweightlist is in EACH session folder.
            try
                alldeweight = [];
                allrepaired = [];
                lastimg = 0; % approach to identify number of scans per session introduced by CP, 2021/05/12; original ArtRepair code assumes equal number of scans per Session
                for isess = 1:num_sess
                    lastimg = lastimg + repSPM.SPM.nscan(isess); % lastimg gets the list number of the last image of this Session. images for all sessions are in  SPM.xY.P, size of session is SPM.nscan(i)   
                    imagedir = fileparts(repSPM.SPM.xY.P(lastimg,:)); % write information on location of volumes to variables % by CP, 2021/05/11
                    deweightlist = fullfile(imagedir,'art_deweighted.txt');
                    repairedlist = fullfile(imagedir,'art_repaired.txt');
                    alldeweight = [alldeweight, load(deweightlist)];
                    allrepaired = [allrepaired, load(repairedlist)];
                end
            catch
                disp('WARNING: Could not locate a deweighting file for each session')
            end
        end
        
        % Collect parameters to characterize improvement
        DeweightedScans_n(subcount,:) = length(alldeweight);
        DeweightedScans_percent(subcount,:) = round (length(alldeweight) / sum(SPM.nscan) * 100 , 2);
        
        RepairedScans_n(subcount,:) = length(allrepaired);
        RepairedScans_percent(subcount,:) = round (length(allrepaired) / sum(SPM.nscan) * 100 ,2);
        
        fprintf(['Finished ',allsubs(sub).name,' successfully!\n'])
        subcount = subcount + 1;
    end
end

% Calculate percent improvement
StdEstimationError_improvement_percent = round ((StdEstimationError_origSPM - StdEstimationError_redoSPM) ./ StdEstimationError_origSPM * 100 , 2);
ResMSerror_improvement_percent = round ((ResMSerror_origSPM - ResMSerror_redoSPM) ./ ResMSerror_origSPM * 100 , 2);

% Round statistics
MeanEstimationError_origSPM = round (MeanEstimationError_origSPM,3);
MeanEstimationError_redoSPM = round (MeanEstimationError_redoSPM,3);
StdEstimationError_origSPM = round(StdEstimationError_origSPM,3);
StdEstimationError_redoSPM = round(StdEstimationError_redoSPM,3);
ResMSerror_origSPM = round(ResMSerror_origSPM,3);
ResMSerror_redoSPM = round(ResMSerror_redoSPM,3);

% Save to summary table
GQsumdir = fullfile(spm_analysis_dir,'firstlevel',sesdir,taskdir,taskspecdir,'GlobalQuality_overview.txt');
datatable = table(datafile_ID,session_data,StdEstimationError_improvement_percent,ResMSerror_improvement_percent,RepairedScans_percent,DeweightedScans_percent,RepairedScans_n,DeweightedScans_n,MeanEstimationError_origSPM,StdEstimationError_origSPM,MeanEstimationError_redoSPM,StdEstimationError_redoSPM,ResMSerror_origSPM,ResMSerror_redoSPM);
writetable(datatable,GQsumdir,'Delimiter','tab');