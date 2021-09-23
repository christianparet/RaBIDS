%% Read data from datasheet
data = readtable('datasheet.xlsx','ReadRowNames',true,'PreserveVariableNames',true,'NumHeaderLines',0);

userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));

DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};

%% provide user input
unzipanat = false; % unzip anatomical niftis? default is false

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
        
        try

            if unzipanat && isfolder([subd,filesep,'anat'])
                zippednii = dir([subd,filesep,'anat',filesep,'sub-*_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz']);

                for anats = 1:length(zippednii)
                    gunzip([subd,filesep,'anat',filesep,zippednii(anats).name])
                    delete([subd,filesep,'anat',filesep,zippednii(anats).name])
                    fprintf(['Unzipped ',zippednii(anats).name,', deleted zipped nii.\n'])
                end
            else
                fprintf('No anat directory found.\n');
            end
            
            allses = dir([subd,filesep,'ses-*']);
            
            for ses = 1:length(allses)
                sesd = [subd,filesep,allses(ses).name];
                fprintf([allses(ses).name,'\n'])
                
                if isfolder([sesd,filesep,'func'])
                    zippednii = dir([sesd,filesep,'func',filesep,'sub-*-preproc_bold.nii.gz']);
                    
                    for funcs = 1:length(zippednii)
                        gunzip([sesd,filesep,'func',filesep,zippednii(funcs).name])
                        delete([sesd,filesep,'func',filesep,zippednii(funcs).name])
                        fprintf(['Unzipped ',zippednii(funcs).name,', deleted zipped nii.\n'])
                    end
                else
                    fprintf('No func directory found.\n');
                end
                    
            end
        catch
            if isfolder([subd,filesep,'func'])
                zippednii = dir([subd,filesep,'func',filesep,'sub-*-preproc_bold.nii.gz']);
                
                for funcs = 1:length(zippednii)
                    gunzip([subd,filesep,'func',filesep,zippednii(funcs).name])
                    delete([sesd,filesep,'func',filesep,zippednii(funcs).name])
                    fprintf(['Unzipped ',zippednii(funcs).name,', deleted zipped nii.\n'])
                end
            else
                fprintf('No func directory found.\n');
            end
        end
    end
end
fprintf('End')