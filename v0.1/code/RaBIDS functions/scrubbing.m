function out = scrubbing(subject,task,data_analysis_path,ses_id,addsub,overwrite)
% v0.1 release

% Scrubbing regressors, adapted from Soroosh Afyouni, NISOx.org, 2017 and adapted from Martin Fungi Gerchens real-time fMRI analysis pipeline
% Framewise displacement (FD) parameters

% %% Comment out if function in use
% clear
% clc
% data_analysis_path = 'S:\AG-Austausch\RaBIDS\example\your project directory';
% subject = 'RaBIDS01';
% ses_id = 'ses-01';
% addsub = 'y'; % add prefix 'sub-' before subject name
% task = 'faces';
% overwrite = 'y';
%

%% Define paths
if contains(addsub,'y')
    prefix = 'sub-';
else
    prefix = '';
end

if strcmp(ses_id,'none')
    write_ses = '_';
    wd = [prefix,subject];
else
    write_ses = ['_',ses_id,'_'];
    wd = [prefix,subject,filesep,ses_id];
end

derivp = fullfile(data_analysis_path,'dataset','derivatives\RaBIDS-prepro\',wd);

rpfile = dir([derivp,filesep,'rp_a',prefix,subject,write_ses,'task-',task,'_bold.txt']);

try
    loc1 = strfind(rpfile.name,'rp_a');
    loc2 = strfind(rpfile.name,'bold.txt');
    nuisancef = [derivp,filesep,'nuisance_',rpfile.name(loc1+4:loc2+3)];
   
catch
    out = 'Realignment parameters not found. Continue with next session.\n';
    return
end

if ~isfile(nuisancef) || contains(overwrite,'y')

    %% Set parameters
    FD_radius = 0.5; % from Fungi; not clear why 0.5 ???
    deg_par = [0 0 0 1 1 1]; % motion parameter columns with degree values
    r_Idx   = find(deg_par);
    FD_threshold = 0.5; % threshold to censor volumes by FD in mm

    %% Search for superthreshold movement and build nuisance matrix
    if ~length(rpfile)==1
        out = 'No or too many realignment parameters found. Continue with next session.\n';
        return
    else
        fprintf(['Processing now ',rpfile.name,'\n'])
        T = readtable([derivp,filesep,rpfile.name]);
        rp = table2array(T);
        MP = rp(:,1:6); % MP = movement parameter
        MP(:,r_Idx) = (2*FD_radius*pi/360)*MP(:,r_Idx);
        dMP  = diff(MP);
        FDts = [0; sum(abs(dMP),2)];
        FD_cens = find(FDts>FD_threshold);

        if any(FD_cens)
            mov_dummy = zeros(size(FDts,1),length(FD_cens));
            for l = 1:length(FD_cens)
                mov_dummy(FD_cens(l),l) = 1;
            end
            R = [rp mov_dummy];
        else
            R = rp;
        end

        save(nuisancef,'R')

        out = 'Scrubbing successful, nuisance regressor written.\n';
    end
    
else
    out = 'Nuisance file exists, continue with next step.\';
end