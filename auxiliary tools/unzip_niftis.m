%% Unzip images in BIDS derivative structure
% Christian Paret, ZI Mannheim, 2021

% Download this program into your dataset/code directory to run

%% Changelog

% 2021/05/20 fixed bug to handle datasets without session-subdirectories
% (i.e., with subject directories having func directories on next lower
% hierarchy.

%% Read data from datasheet
data = readtable('datasheet.xlsx','ReadRowNames',true,'PreserveVariableNames',true,'NumHeaderLines',0);

userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));

DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};

%% unzip images
derivd = dir([data_analysis_path,filesep,'derivatives',filesep,'fmriprep*']);

fprintf(['Found ',num2str(length(derivd)),' derivative directories:\n'])

for i=1:length(derivd)
    fprintf([num2str(i),') ',derivd(i).name,'\n'])
end

selectd = input('Enter number of derivative directory to unzip preporcessed niftis and press enter.\n');

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

allsubs = dir([derivd,filesep,'sub-*']);

for sub = 1:length(allsubs)
    subd = [derivd,filesep,allsubs(sub).name];
    
    if isfolder(subd)
        
        fprintf([allsubs(sub).name,'\n'])
        
        allses = dir([subd,filesep,'ses-*']);        
        nses = length(allses);
        if nses < 1 % Account for datasets that do not have ses subdirectory
            nses = 1;
            sesid = false;
        else
            sesid = true;
        end
        
        for ses = 1:nses
            if sesid
                sesdir = [subd,filesep,allses(ses).name];
                fprintf([allses(ses).name,'\n'])
            else
                sesdir = subd;
            end
            
            if isfolder([sesdir,filesep,'func'])
                zippednii = dir([sesdir,filesep,'func',filesep,'sub-*-preproc_bold.nii.gz']);
                
                if isempty(zippednii)
                    fprintf('No zipped images found.\n')
                else
                    for funcs = 1:length(zippednii)
                        fprintf(['Unzipping ',zippednii(funcs).name,'...\n'])
                        gunzip([sesdir,filesep,'func',filesep,zippednii(funcs).name])
                        fprintf('Unzip successful, deleting zipped nii...\n')
                        delete([sesdir,filesep,'func',filesep,zippednii(funcs).name])
                        fprintf('Success!\n')
                    end
                end
            else
                fprintf('No func directory found.\n')
            end
        end
    end
end
fprintf('End')