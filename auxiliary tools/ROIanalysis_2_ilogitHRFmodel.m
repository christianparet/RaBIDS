%% ROIanalysis_2_HRFmodel.m
% Christian Paret, ZI Mannheim, 2022
%
% Estimate HRF parameters using the inverse logit model (Lindquist & Wager, 2007)
% Ref: Lindquist, M.A., Wager, T.D., 2007. Validity and Power in Hemodynamic Response Modeling: A Comparison Study and a New Approach. Hum Brain Mapp 28, 764–784. https://doi.org/10.1002/hbm.20310
% Below program is based on template code which comes with the CanlabCore package (https://canlab.github.io/). The template code can be accessed here: CanlabCore\HRF_Est_Toolbox2\Example.m
%
% Run this program on BOLD timecourse data
% To work properly you should have your data arranged with RaBIDS (https://github.com/christianparet/RaBIDS).
% First-level data should be prepared with the program SPManalysis_1_firstlevel.m, and BOLD signal data should be extracted with ROIanalysis_1_extract_eigenvariate.m. Both programs are available in the RaBIDS auxiliary tools.
% Program expects to find BOLD signal data for the first subject in a file structure like this: E:\mytrainingdata\your project directory\spm_analysis\firstlevel\ses-01\task-faces\taskrelated_activity\original\sub-RABIDS01_ses-01_task-faces.
% The matlab structure containing BOLD signal data is named sth like 'VOI_mask-RightAmygdala25_sess_1.mat'
%
% Cave: Program assumes that stimulus onsets are defined in unit seconds.
%       Status 2022/05/16: Double-check of upsampling approach remains pending
%       Check whether different resampleTRs as well as interpolation methods result in consistent HRF parameter outcomes

%%
clear
clc
diary off

%% Add path to CanlabCore
addpath(genpath('E:\CanlabCore\CanlabCore'));

%% Settings
VOIname = 'Right Amygdala 25'; % enter VOI name; e.g. if file name is VOI_mask-visualcortex_sess_1.mat: VOIname = 'visualcortex';

origTR = 2; % original TR in seconds
resampleTR = 0.1; % resample to TR (seconds), should be in the order of stimulus onset time (SOT) resolution
interpolmeth = 'spline'; % spline, pchip or makima. Not clear which one is better
getFIR = 1; % FIR fitting fails with resampleTR<0.01
T = 30; % unit=seconds
FWHM = 4; % FWHM for residual scan
pval = 0.01;
df = 600;
alpha = 0.001;

%% Build HRF
[xBF] = spm_get_bf(struct('dt', resampleTR, 'name', 'hrf (with time and dispersion derivatives)', 'length', 32));
clear Xtrue
for i = 1:1, xx = conv(xBF.bf(:,i), [1 1 1 1 1 1 ]');
    Xtrue(:, i) = xx(1:32/resampleTR);
end
for i = 2:3, xx = conv(xBF.bf(:,i), [1 1]');
    Xtrue(:, i) = xx(1:32/resampleTR);
end
hrf = Xtrue * [1 .3 .2]';
xsecs = 0:resampleTR:32;

hrf = [ 0; 0; hrf];
hrf = hrf(1:length(xsecs));
hrf = hrf ./ max(hrf);

%% Define directories and initiate variables
fprintf('Select directory containing subject directories with SPM.mat files.\nSomething like %s.','\mytrainingdata\your project directory\spm_analysis\firstlevel\ses-01\task-faces\taskrelated_activity\original')
firstleveldir = uigetdir(pwd,'Select firstlevel directory that includes the subject directories. See command window for help.'); % directory containing child directories with each child containing an estimated firstlevel SPM.mat
% firstleveldir = 'Y:\Projects\EFPTest\Data_analysis\spm_analysis\firstlevel\ses-pre\task-efptest\taskrelated_activity\repaired';
subdirs = dir([firstleveldir,filesep,'sub*']);

diaryfn = 'diary_ROIanalysis_2_ilogitHRFmodel.txt';
diary(fullfile(firstleveldir,diaryfn)) % Write command window output to diary logfile

subject_count = 1;

%% Write some output to command window
diary on
fprintf(['\n\n--------------- ROI analysis: Estimate ilogit HRF model (Lindquist & Wager, 2007) ---------------\n\nExtracted BOLD signal timecourse is from volume of interest: ',VOIname,'\nNote: Program can only handle SPM models with a single scanning session.\n\n']);

%%
if ~isempty(subdirs)
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

            % Locate timecourse files
            tcfile = fullfile(firstleveldir,subdirs(sub).name,['VOI_mask-',VOIname,'_sess_1.mat']);

            % Locate firstlevel spm.mat
            spm_name = fullfile(firstleveldir,subdirs(sub).name,'SPM.mat');

            if exist(spm_name,'file') && exist(tcfile,'file')
                
                % write SubjectID vector for table output of group results
                SubjectID{subject_count,:} = subdirs(sub).name(1:dum1-1); 

                % Load BOLD timecourse to workspace
                VOIdata = load(tcfile); % load timecourse to variable VOIdata

                %% Resample BOLD timecourse to target resolution
                % Interpolate BOLD timecourse and resample
                x = 1:1:length(VOIdata.Y);
                xx = 1:resampleTR/origTR:length(VOIdata.Y); % determine new resolution of x-axis
                switch interpolmeth
                    case 'spline'
                        yy = spline(x,VOIdata.Y,xx);
                    case 'pchip'
                        yy = pchip(x,VOIdata.Y,xx);
                    case 'makima'
                        yy = makima(x,VOIdata.Y,xx);
                    otherwise
                        fprintf('Signal interpolation method is invalid. Check value of variable resmeth.\n')
                        break
                end


%                 plot(x,VOIdata.Y,'o',xx,yy); % Use this command to visually inspect the interpolation results

                % Get the stimulus onset times
                spmfile = load(spm_name);

                %% Make output directory to save results to
                outdir = fullfile(firstleveldir,subdirs(sub).name,VOIname,['InterpolateMethod_',interpolmeth]); % directory to save results to

                %% Check for existing file with estimated parameters
                usepa = 1; % variable usepa controls for existing parameter data
                if exist(fullfile(outdir,'HRFmodeling.mat')) % was the model estimated before?
                    fprintf('    Found existing parameters logfile ''HRFmodeling.mat''.\n')
                    prelim = load(fullfile(outdir,'HRFmodeling.mat')); % load model and check consistency with parameters
%                     if ~strcmp(prelim.parameters.subID,subdirs(sub).name(1:dum1-1))
%                         fprintf('     --> wrong subject ID, overwriting existing parameters logfile.\n')
%                         usepa = 'false';
%                     end
                    if ~strcmp(prelim.parameters.VOIname,VOIname)
                        fprintf('     --> different VOI name, overwriting existing parameters logfile.\n')
                        usepa = 0;
                    end
                    if prelim.parameters.OriginalTR~=origTR
                        fprintf('     --> different original TR, overwriting existing parameters logfile.\n')
                        usepa = 0;
                    end
                    if prelim.parameters.InterpolateToTR~=resampleTR
                        fprintf('     --> interpolated to different TR, overwriting existing parameters logfile.\n')
                        usepa = 0;
                    end
                    if prelim.parameters.HRFLength~=T
                        fprintf('     --> different HRF length ''T'', overwriting existing parameters logfile.\n')
                        usepa = 0;
                    end
                    if prelim.parameters.FWHM~=FWHM
                        fprintf('     --> different FWHM, overwriting existing parameters logfile.\n')
                        usepa = 0;
                    end
                    if length(spmfile.SPM.Sess.U)~=length(prelim.parameters.conditions)
                        fprintf('     --> different number of conditions, overwriting existing parameters logfile.\n')
                        usepa = 0;
                    else
                        for cond = 1:length(spmfile.SPM.Sess.U)
                            if ~strcmp(prelim.parameters.conditions{cond},spmfile.SPM.Sess.U(cond).name)
                                fprintf('     --> different condition name, overwriting existing parameters logfile.\n')
                                usepa = 0;
                                break
                            end
                        end
                    end
                    
                else
                    mkdir(outdir); % assuming directory is missing
                    usepa = 0;
                end

                %% Retrieve or estimate HRF parameters
                clear parameters % parameters from HRF modeling have been/will be written to structure called 'parameters'

                if usepa
                    fprintf('    Existing meta-data is consistent with current settings.\n    Retrieving parameters....\n')

                    % Feed matlab-structure "parameters" with existing subject data
                    parameters.subID = subdirs(sub).name(1:dum1-1); parameters.VOIname=prelim.parameters.VOIname; 
                    parameters.OriginalTR=prelim.parameters.OriginalTR;
                    parameters.InterpolateToTR=prelim.parameters.InterpolateToTR;
                    parameters.InterpolateMethod=prelim.parameters.InterpolateMethod; 
                    parameters.HRFLength=prelim.parameters.HRFLength;
                    parameters.FWHM=prelim.parameters.FWHM;
                    parameters.conditions=prelim.parameters.conditions;
                    parameters.ilogitHRF_AUC=prelim.parameters.ilogitHRF_AUC; 
                    parameters.ilogitHRF_estimatedHRF=prelim.parameters.ilogitHRF_estimatedHRF; 
                    parameters.ilogitHRF_Amplitude=prelim.parameters.ilogitHRF_Amplitude;
                    parameters.ilogitHRF_TimeToPeak=prelim.parameters.ilogitHRF_TimeToPeak;
                    parameters.ilogitHRF_Width=prelim.parameters.ilogitHRF_Width;
                    parameters.ilogitHRF_Mismodeling=prelim.parameters.ilogitHRF_Mismodeling; 
                    parameters.canonicalHRFandDV_AUC=prelim.parameters.canonicalHRFandDV_AUC; 
                    parameters.canonicalHRFandDV_estimatedHRF=prelim.parameters.canonicalHRFandDV_estimatedHRF;
                    parameters.canonicalHRFandDV_Amplitude=prelim.parameters.canonicalHRFandDV_Amplitude; 
                    parameters.canonicalHRFandDV_TimeToPeak=prelim.parameters.canonicalHRFandDV_TimeToPeak;
                    parameters.canonicalHRFandDV_Width=prelim.parameters.canonicalHRFandDV_Width; 
                    parameters.canonicalHRFandDV_Mismodeling=prelim.parameters.canonicalHRFandDV_Mismodeling;

                    if getFIR
                        parameters.sFIR_AUC=prelim.parameters.sFIR_AUC; 
                        parameters.sFIR_estimatedHRF=prelim.parameters.sFIR_estimatedHRF;
                        parameters.sFIR_Amplitude=prelim.parameters.sFIR_Amplitude; 
                        parameters.sFIR_TimeToPeak=prelim.parameters.sFIR_TimeToPeak; 
                        parameters.sFIR_Width=prelim.parameters.sFIR_Width;
                        parameters.sFIR_Mismodeling=prelim.parameters.sFIR_Mismodeling; 
                    end

                else
                    % If parameters have not been estimated yet: estimate HRF
                    fprintf('    Results of HRF modeling will be saved to: %s',outdir), fprintf('\n')

                    % Feed matlab-structure "parameters" with meta-data defined above for HRF estimation
                    parameters.subID = subdirs(sub).name(1:dum1-1); parameters.VOIname=VOIname; parameters.OriginalTR=origTR; parameters.InterpolateToTR=resampleTR; parameters.InterpolateMethod=interpolmeth; parameters.HRFLength=T; parameters.FWHM=FWHM;

                    for cond = 1:length(spmfile.SPM.Sess.U) % run through all conditions
                        parameters.conditions(cond) = spmfile.SPM.Sess.U(cond).name;
                        fprintf('    Processing condition ''%s''',parameters.conditions{cond}), fprintf('\n')

                        % Visualize result of interpolation
                        [f1,~] = create_figure; subplot(3,1,1);
                        % han = plot(x,VOIdata.Y,'o',xx,yy); % to display original signal on resampled signal
                        han = plot(yy);
                        title([subdirs(sub).name,' ',VOIname]); drawnow
                        hold on;
                        Run = zeros(length(yy),1);
                        resampleons = round( spmfile.SPM.Sess.U(cond).ons / resampleTR ); % change resolution of onsets to be consistent with resampled BOLD timecourse
                        if resampleons(1)==0, resampleons(1)=1; end % zero-entries are invalid
                        for i=1:length(spmfile.SPM.Sess.U(cond).ons), Run(resampleons(i)) = 1; end
                        hh = plot_onsets(Run,'k',-3,1, 1/resampleTR);
                        drawnow
                        Runc{1} = Run;

                        %% Fit HRF using IL-function
                        % Choose mode (deterministic/stochastic)
                        mode = 0;   % 0 - deterministic aproach
                        % 1 - simulated annealing approach
                        % Please note that when using simulated annealing approach you
                        % may need to perform some tuning before use.

                        [h1, fit1, e1, param] = Fit_Logit2(yy',resampleTR,Runc,T,mode);
                        [pv sres sres_ns1] = ResidScan(e1, FWHM);
                        [PowLoss1] = PowerLoss(e1, fit1, (length(yy)-7) , yy', resampleTR, Runc, alpha);

                        hold on; han(2) = plot(fit1,'r');

                        disp('Summary: IL_function');

                        disp('Amplitude:'); disp(param(1));
                        disp('Time-to-peak:'); disp(param(2)*resampleTR);
                        disp('Width:'); disp(param(3)*resampleTR);

                        % calculate Area Under the Curve (AUC)
                        parameters.ilogitHRF_AUC(cond) = trapz(h1)*resampleTR; % Codes needs to be verified!
                        parameters.ilogitHRF_estimatedHRF(cond,:) = h1;
                        parameters.ilogitHRF_Amplitude(cond) = param(1); parameters.ilogitHRF_TimeToPeak(cond) = param(2)*resampleTR; parameters.ilogitHRF_Width(cond) = param(3)*resampleTR; parameters.ilogitHRF_Mismodeling(cond) = pv;

                        disp('MSE:'); disp((1/(length(yy)-1)*sum(e1.^2)));
                        disp('Mis-modeling:'); disp(pv);
                        disp('Power Loss:'); disp(PowLoss1);

                        %% Fit HRF using FIR-model
                        hancount = 3;
                        if getFIR
                            % Choose mode (FIR/sFIR)
                            mode = 1;   % 0 - FIR
                            % 1 - smooth FIR

                            [h2, fit2, e2, param] = Fit_sFIR(yy',resampleTR,Runc,T,mode);
                            [pv sres sres_ns2] = ResidScan(e2, FWHM);
                            [PowLoss2] = PowerLoss(e2, fit2, (length(yy)-T) , yy', resampleTR, Runc, alpha);

                            hold on; han(hancount) = plot(fit2,'g');

                            disp('Summary: FIR');

                            disp('Amplitude'); disp(param(1));
                            disp('Time-to-peak'); disp(param(2)*resampleTR);
                            disp('Width'); disp(param(3)*resampleTR);

                            % calculate Area Under the Curve (AUC)
                            parameters.sFIR_AUC(cond) = trapz(h2)*resampleTR; % Codes needs to be verified!
                            parameters.sFIR_estimatedHRF(cond,:) = h2;
                            parameters.sFIR_Amplitude(cond) = param(1); parameters.sFIR_TimeToPeak(cond) = param(2)*resampleTR; parameters.sFIR_Width(cond) = param(3)*resampleTR; parameters.sFIR_Mismodeling(cond) = pv;

                            disp('MSE:'); disp((1/(length(yy)-1)*sum(e2.^2)));
                            disp('Mis-modeling'); disp(pv);
                            disp('Power Loss:'); disp(PowLoss2);

                            labels = {'Data' 'IL' 'sFIR' 'DD'};
                            hancount=hancount+1;
                        else
                            labels = {'Data' 'IL' 'DD'};
                        end

                        %% Fit HRF using Canonical HRF + 2 derivatives
                        p=1;
                        [h3, fit3, e3, param, out] = Fit_Canonical_HRF(yy',resampleTR,Runc,30,p);
                        [pv sres sres_ns3] = ResidScan(e3, FWHM);
                        [PowLoss3] = PowerLoss(e3, fit3, (length(yy)-p) , yy', resampleTR, Runc, alpha);

                        hold on; han(hancount) = plot(fit3,'m');

                        legend(han,labels)

                        disp('Summary: Canonical + 2 derivatives');

                        disp('Amplitude'); disp(param(1));
                        disp('Time-to-peak'); disp(param(2)*resampleTR);
                        disp('Width'); disp(param(3)*resampleTR);

                        % calculate Area Under the Curve (AUC)
                        parameters.canonicalHRFandDV_AUC(cond) = trapz(h3)*resampleTR; % Codes needs to be verified!
                        parameters.canonicalHRFandDV_estimatedHRF(cond,:) = h3;
                        parameters.canonicalHRFandDV_Amplitude(cond) = param(1); parameters.canonicalHRFandDV_TimeToPeak(cond) = param(2)*resampleTR; parameters.canonicalHRFandDV_Width(cond) = param(3)*resampleTR; parameters.canonicalHRFandDV_Mismodeling(cond) = pv;

                        disp('MSE:'); disp((1/(length(yy)-1)*sum(e3.^2)));
                        disp('Mis-modeling'); disp(pv);
                        disp('Power Loss:'); disp(PowLoss3);

                        %% Plot all of this
                        subplot(3,2,5); hold on;
                        plot(xsecs, hrf, 'k')
                        xsecs1 = xsecs(1:length(h1));
                        han2 = plot(xsecs1, h1,'r');
                        dum = 2;
                        if getFIR
                            xsecs2 = xsecs(1:length(h2));
                            han2(dum) = plot(xsecs2, h2,'g');
                            labels = {'IL' 'sFIR' 'DD'};
                            dum=dum+1;
                        else
                            labels = {'IL' 'DD'};
                        end
                        xsecs3 = xsecs(1:length(h3));
                        han2(dum) = plot(xsecs3, h3,'m');
                        legend(han2,labels)
                        title('Estimated HRF');

                        subplot(3,1,2); hold on;
                        hh; drawnow

                        han3 = plot(sres_ns1,'r');
                        dum=2;
                        if getFIR
                            hold on; han3(dum) = plot(sres_ns2,'g');
                            labels = {'IL' 'sFIR' 'DD'};
                            dum=dum+1;
                        else
                            labels = {'IL' 'DD'};
                        end
                        hold on; han3(dum) = plot(sres_ns3,'m');
                        hold on; plot((1:length(yy)),zeros(length(yy),1),'--k');
                        legend(han3,labels)
                        title('Mis-modeling (time course)');

                        subplot(3,2,6); hold on;

                        [s1] = Fit_sFIR(sres_ns1,resampleTR,Runc,T,0);
                        if getFIR
                            [s2] = Fit_sFIR(sres_ns2,resampleTR,Runc,T,0);
                        end
                        [s3] = Fit_sFIR(sres_ns3,resampleTR,Runc,T,0);

                        han4 = plot(s1(1:T),'r');
                        dum=2;
                        if getFIR
                            hold on; han4(dum) = plot(s2(1:T),'g');
                            labels = {'IL' 'sFIR' 'DD'};
                            dum=dum+1;
                        else
                            labels = {'IL' 'DD'};
                        end
                        hold on; han4(dum) = plot(s3(1:T),'m');
                        hold on; plot((1:T),zeros(T,1),'--k');
                        legend(han4,labels)
                        title('Mis-modeling (HRF)');

                        hold off

                        %% Save results
                        savefig(f1,fullfile(outdir,['HRFmodeling_',spmfile.SPM.Sess.U(cond).name{:},'.fig']))
                    end
                    save(fullfile(outdir,'HRFmodeling.mat'),'parameters')
                end

                %% Writing parameter output table
                % Run through all conditions in order to create a string variable listing all parameter variables. Create string-variable which is an executable command for producing a vector with all parameters included. Use eval command to execute the string as a matlab command.
                for cond = 1:length(parameters.conditions)
                    % Inverse logit HRF parameters
                    eval(sprintf('%s_IL_Amplitude(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.ilogitHRF_Amplitude(cond)));
                    eval(sprintf('%s_IL_TimeToPeak(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.ilogitHRF_TimeToPeak(cond)));
                    eval(sprintf('%s_IL_Width(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.ilogitHRF_Width(cond)));
                    eval(sprintf('%s_IL_AUC(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.ilogitHRF_AUC(cond)));
                    eval(sprintf('%s_IL_Mismodeling(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.ilogitHRF_Mismodeling(cond)));

                    % Canonical HRF + 2 derivative parameters
                    eval(sprintf('%s_canonical_Amplitude(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.canonicalHRFandDV_Amplitude(cond)));
                    eval(sprintf('%s_canonical_TimeToPeak(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.canonicalHRFandDV_TimeToPeak(cond)));
                    eval(sprintf('%s_canonical_Width(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.canonicalHRFandDV_Width(cond)));
                    eval(sprintf('%s_canonical_AUC(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.canonicalHRFandDV_AUC(cond)));
                    eval(sprintf('%s_canonical_Mismodeling(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.canonicalHRFandDV_Mismodeling(cond)));

                    if getFIR
                        % sFIR HRF parameters
                        eval(sprintf('%s_sFIR_Amplitude(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.sFIR_Amplitude(cond)));
                        eval(sprintf('%s_sFIR_TimeToPeak(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.sFIR_TimeToPeak(cond)));
                        eval(sprintf('%s_sFIR_Width(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.sFIR_Width(cond)));
                        eval(sprintf('%s_sFIR_AUC(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.sFIR_AUC(cond)));
                        eval(sprintf('%s_sFIR_Mismodeling(subject_count,1) = %.4f;', parameters.conditions{cond},parameters.sFIR_Mismodeling(cond)));
                    end

%                     eval([parameters.conditions{cond},'_IL_Amplitude(subject_count,1) = ', num2str(parameters.ilogitHRF_Amplitude(cond)),';']);
%                     eval([parameters.conditions{cond},'_IL_TimeToPeak(subject_count,1) = ', num2str(parameters.ilogitHRF_TimeToPeak(cond)),';']);
%                     eval([parameters.conditions{cond},'_IL_Width(subject_count,1) = ', num2str(parameters.ilogitHRF_Width(cond)),';']);
%                     eval([parameters.conditions{cond},'_IL_AUC(subject_count,1) = ', num2str(parameters.ilogitHRF_AUC(cond)),';']);
%                     eval([parameters.conditions{cond},'_IL_Mismodeling(subject_count,1) = ', num2str(parameters.ilogitHRF_Mismodeling(cond)),';']);
%                     
%                     eval([parameters.conditions{cond},'_canonical_Amplitude(subject_count,1) = ', num2str(parameters.canonicalHRFandDV_Amplitude(cond)),';']);
%                     eval([parameters.conditions{cond},'_canonical_TimeToPeak(subject_count,1) = ', num2str(parameters.canonicalHRFandDV_TimeToPeak(cond)),';']);
%                     eval([parameters.conditions{cond},'_canonical_Width(subject_count,1) = ', num2str(parameters.canonicalHRFandDV_Width(cond)),';']);
%                     eval([parameters.conditions{cond},'_canonical_AUC(subject_count,1) = ', num2str(parameters.canonicalHRFandDV_AUC(cond)),';']);
%                     eval([parameters.conditions{cond},'_canonical_Mismodeling(subject_count,1) = ', num2str(parameters.canonicalHRFandDV_Mismodeling(cond)),';']);

%                     if getFIR
%                         % sFIR HRF parameters
%                         eval([parameters.conditions{cond},'_sFIR_Amplitude(subject_count,1) = ', num2str(parameters.sFIR_Amplitude(cond)),';']);
%                         eval([parameters.conditions{cond},'_sFIR_TimeToPeak(subject_count,1) = ', num2str(parameters.sFIR_TimeToPeak(cond)),';']);
%                         eval([parameters.conditions{cond},'_sFIR_Width(subject_count,1) = ', num2str(parameters.sFIR_Width(cond)),';']);
%                         eval([parameters.conditions{cond},'_sFIR_AUC(subject_count,1) = ', num2str(parameters.sFIR_AUC(cond)),';']);
%                         eval([parameters.conditions{cond},'_sFIR_Mismodeling(subject_count,1) = ', num2str(parameters.sFIR_Mismodeling(cond)),';']);
%                     end

                    % Write estimated HRF timecourse: write each time bin to separate varialble
                    for bin = 1:length(parameters.ilogitHRF_estimatedHRF) % assuming all HRFs have same length
                        eval(sprintf('%s_IL_bin_%04d (subject_count,1) = %.4f;',parameters.conditions{cond},bin,parameters.ilogitHRF_estimatedHRF(cond,bin)));
                        eval(sprintf('%s_canonical_bin_%04d (subject_count,1) = %.4f;',parameters.conditions{cond},bin,parameters.canonicalHRFandDV_estimatedHRF(cond,bin)));
                        if getFIR
                            eval(sprintf('%s_sFIR_bin_%04d (subject_count,1) = %.4f;',parameters.conditions{cond},bin,parameters.sFIR_estimatedHRF(cond,bin)));
                        end
                    end
                end
                subject_count = subject_count + 1;
            else
                fprintf('    No SPM.mat and/or no timecourse file (VOI_mask-*.mat) found for this subject.\n')
            end
        end
    end

    %% Write all subjects to output table
    getparamtable = 'SubjectID';
    getILtable = 'SubjectID';
    getcanonicaltable = 'SubjectID';

    outdir2 = fullfile(firstleveldir,'estimated_HRF',VOIname,['InterpolateMethod_',interpolmeth]);
    fprintf('    Saving group results to %s\n',outdir2)

    if exist(outdir2,'dir')
        fprintf('    Output directory for group results already exists.\n    Overwriting existing files...\n')
    else
        mkdir(outdir2)
    end

    fprintf(['    Command window output is written to file ',diaryfn,' and is saved next to the subject directories.\n'])

    if getFIR
        getsFIRtable = 'SubjectID';
        % Run through all conditions in order to create a string variable listing all parameter variables
        for cond = 1:length(parameters.conditions)
            % Create string-variable which is an executable command for producing a vector with all parameters included. Use eval command to execute the string as a matlab command.
            eval(sprintf('getparamtable = [getparamtable,'', %s_IL_Amplitude'','', %s_IL_TimeToPeak'','', %s_IL_Width'','', %s_IL_AUC'','', %s_IL_Mismodeling'','', %s_canonical_Amplitude'','', %s_canonical_TimeToPeak'','', %s_canonical_Width'','', %s_canonical_AUC'','', %s_canonical_Mismodeling'','', %s_sFIR_Amplitude'','', %s_sFIR_TimeToPeak'','', %s_sFIR_Width'','', %s_sFIR_AUC'','', %s_sFIR_Mismodeling''];',...
                parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond}));
            for bin = 1:length(parameters.ilogitHRF_estimatedHRF) % assuming all HRFs have same length
                % Create string-variable which is an executable command for producing a vector with all parameters included. Use eval command to execute the string as a matlab command.
                eval(sprintf('getILtable = [getILtable,'', %s_IL_bin_%04d''];',parameters.conditions{cond},bin));
                eval(sprintf('getcanonicaltable = [getcanonicaltable,'', %s_canonical_bin_%04d''];',parameters.conditions{cond},bin));
                eval(sprintf('getsFIRtable = [getsFIRtable,'', %s_sFIR_bin_%04d''];',parameters.conditions{cond},bin));
            end
        end
        eval(['sFIR_estimated_tc = table(',getsFIRtable,');']);
        writetable(sFIR_estimated_tc,fullfile(outdir2,'sFIR_estimated_tc.txt'),'Delimiter','tab');
    else
        % Run through all conditions in order to create a string variable listing all parameter variables
        for cond = 1:length(parameters.conditions)
            % Create string-variable which is an executable command for producing a vector with all parameters included. Use eval command to execute the string as a matlab command.
            eval(sprintf('getparamtable = [getparamtable,'', %s_IL_Amplitude'','', %s_IL_TimeToPeak'','', %s_IL_Width'','', %s_IL_AUC'','', %s_IL_Mismodeling'','', %s_canonical_Amplitude'','', %s_canonical_TimeToPeak'','', %s_canonical_Width'','', %s_canonical_AUC'','', %s_canonical_Mismodeling''];',...
                parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond},parameters.conditions{cond}));
            for bin = 1:length(parameters.ilogitHRF_estimatedHRF) % assuming all HRFs have same length
                % Create string-variable which is an executable command for producing a vector with all parameters included. Use eval command to execute the string as a matlab command.
                eval(sprintf('getILtable = [getILtable,'', %s_IL_bin_%04d''];',parameters.conditions{cond},bin));
                eval(sprintf('getcanonicaltable = [getcanonicaltable,'', %s_canonical_bin_%04d''];',parameters.conditions{cond},bin));
            end
        end

    end
    % Use table command to create a table of parameter vectors
    eval(['HRF_parameters = table(',getparamtable,');']);
    eval(['IL_estimated_tc = table(',getILtable,');']);
    eval(['Canonical_estimated_tc = table(',getcanonicaltable,');']);

    % Save table as text files to output directory
    writetable(HRF_parameters,fullfile(outdir2,'HRF_parameters.txt'),'Delimiter','tab');
    writetable(IL_estimated_tc,fullfile(outdir2,'IL_estimated_tc.txt'),'Delimiter','tab');
    writetable(Canonical_estimated_tc,fullfile(outdir2,'Canonical_estimated_tc.txt'),'Delimiter','tab');

    % Write a readme file to the output directory
    title = 'Data in HRF parameter-files and estimated tc-files produced with HRF_Est_Toolbox2';
    ref_lit = 'Ref: Lindquist, M.A., Wager, T.D., 2007. Validity and Power in Hemodynamic Response Modeling: A Comparison Study and a New Approach. Hum Brain Mapp 28, 764–784. https://doi.org/10.1002/hbm.20310';
    ref_web = 'Code is based on file CanlabCore\HRF_Est_Toolbox2\Example.m provided under https://canlab.github.io/';
    met1 = 'IL: Parameters received using the inverse-logit model by Lindquist & Wager (2007). Canonical: Canonical HRF model from SPM + 2 derivatives were used.\nsFIR: smoothed FIR model.';
    met2 = ['Original TR of fMRI data: ',num2str(origTR),'. Upsampled to: ',num2str(resampleTR),'. Interpolation method for upsampling: ',interpolmeth,'.'];

    filename = fullfile(outdir2,'readme.txt');
    readme_file = fopen(filename,'w');
    fprintf(readme_file,'%s\n\n%s\n%s\n%s\n\n%s\n%s', title, ref_web, met1, met2,'References:', ref_lit);
    fclose(readme_file);

else
    fprintf('    No subject directories found in selected directory.\n')
end

fprintf('End.\n')
diary off