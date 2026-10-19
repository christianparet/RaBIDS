%% Read data from datasheet
% Program to unzip nifitis following preprocessing with fMRIPrep. Unzipped
% niftis are required for standard SPM firstlevel analysis.
% Prerequisite is the existence of a RaBIDS datasheet excel-file.

% Program unzips functional data. If you want to unzip anatomical data,
% too, you need to set variable 'unzipanat' below to 'true'.
%%

data = readtable('datasheet.xlsx','ReadRowNames',true,'PreserveVariableNames',true,'NumHeaderLines',0);

userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));

DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};

%% provide user input
unzipanat = false; % unzip anatomical niftis? default is false

%% Identify derivative directory to work on
derivd = uigetdir(pwd,'Select subject directory with derivatives (''e.g.: //derivatives/sub-XXX''.');

allsubs = dir([derivd,filesep,'sub-*']);
nsub = length(allsubs);
for sub = 1:nsub

    allses = dir(fullfile(derivd,allsubs(sub).name,'ses-*'));
    nses = length(allses);
    if nses < 1 % Account for datasets that do not have ses subdirectory
        nses = 1;
        sesid = false;
    else
        sesid = true;
    end
    
    for ses = 1:nses
        if sesid
            sesdir = fullfile(derivd,allsubs(sub).name,allses(ses).name);
        else
            sesdir = fullfile(derivd,allsubs(sub).name,allsubs(sub).name);
        end
    
        try
            if unzipanat && isfolder([sesdir,filesep,'anat'])
                zippednii = dir([sesdir,filesep,'anat',filesep,'sub-*_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz']);
    
                for anats = 1:length(zippednii)
                    fprintf('Processing...\n')
                    gunzip([sesdir,filesep,'anat',filesep,zippednii(anats).name])
                    delete([sesdir,filesep,'anat',filesep,zippednii(anats).name])
                    fprintf(['Unzipped ',zippednii(anats).name,', deleted zipped nii.\n'])
                end
            end
    
            if isfolder([sesdir,filesep,'func'])
                zippednii = dir([sesdir,filesep,'func',filesep,'sub-*-preproc_bold.nii.gz']);
    
                for funcs = 1:length(zippednii)
                    fprintf('Processing...\n')
                    gunzip([sesdir,filesep,'func',filesep,zippednii(funcs).name])
                    delete([sesdir,filesep,'func',filesep,zippednii(funcs).name])
                    fprintf(['Unzipped ',zippednii(funcs).name,', deleted zipped nii.\n'])
                end
            else
                fprintf('No func directory found.\n');
            end
    
        catch
            fprintf('Some problem occured./')
        end
    
    end

end

