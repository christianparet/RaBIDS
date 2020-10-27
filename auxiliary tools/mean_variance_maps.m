%% Produce mean and variance images from BOLD series
% Computes mean and variance image from MR series in nifti format and saves it to the data analysis path. This is good for MRI quality control of T2* sequences for MR BOLD imaging. For more information on QC see https://imaging.mrc-cbu.cam.ac.uk/imaging/DataDiagnostics
% SPM path needs to be set. Tested with SPM12 and Matlab v2020. Assumes dataset in BIDS format.
% Christian Paret, ZI Mannheim, 2020

clear
% spm_jobman('initcfg')

%% Select dataset
% Use this routine if dataset has been created with RaBIDS

data = readtable('datasheet.xlsx','ReadRowNames',true,'PreserveVariableNames',true,'NumHeaderLines',0);

userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));

DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};

datasetd = [data_analysis_path,filesep,'dataset'];

% Otherwise: define parent directory of dataset in variable 'data_analysis_path' 

%% 

allsubs = dir([datasetd,filesep,'sub-*']);

for sub = 1:length(allsubs)
    subd = [datasetd,filesep,allsubs(sub).name];
    
    if isfolder(subd)
        
        fprintf([allsubs(sub).name,'\n'])
        
        allses = dir([subd,filesep,'ses-*']);
        
        if ~isempty(allses)
            nses = length(allses);
        else
            allses(1).name = ['.',filesep];
            nses = 1;
        end
        
        for ses = 1:nses
            sesd = [subd,filesep,allses(ses).name];
            if nses == 1
                fprintf('One scan session.\n')
            else
                fprintf([allses(ses).name,'\n'])
            end
            
            if isfolder([sesd,filesep,'func'])
                
                boldnii = dir([sesd,filesep,'func',filesep,'sub-*_bold.nii']);
                boldniigz = dir([sesd,filesep,'func',filesep,'sub-*_bold.nii.gz']);
                
                if ~isempty(boldnii) || ~isempty(boldniigz)
                    if ~isfolder(fullfile(data_analysis_path,'QC',allsubs(sub).name))
                        mkdir(fullfile(data_analysis_path,'QC',allsubs(sub).name))
                    end
                end
                
                if ~isempty(boldnii) && ~isempty(boldniigz)
                    go = input('Found both zipped and unzipped nifti images.\nUnzip may overwrite existing images. Continue? [''Y''/''N'']\n');
                    if strcmp(go,'N')
                        return
                    end
                end
                
                for funcs = 1:length(boldniigz)
                    [~,boldniiunzip(funcs).name,~] = fileparts(boldniigz(funcs).name);
                    gunzip([sesd,filesep,'func',filesep,boldniigz(funcs).name])
                end
                
                boldnii = dir([sesd,filesep,'func',filesep,'sub-*_bold.nii']);
                
                for funcs = 1:length(boldnii)
                    
                    matlabbatch{1}.spm.util.imcalc.input = cellstr(spm_select('expand',fullfile(sesd,'func',boldnii(funcs).name)));
                    matlabbatch{1}.spm.util.imcalc.output = [boldnii(funcs).name(1:end-4),'_mean.nii'];
                    matlabbatch{1}.spm.util.imcalc.outdir = {fullfile(data_analysis_path,'QC',allsubs(sub).name)};
                    matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)';
                    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
                    matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
                    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
                    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
                    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
                    
                    matlabbatch{2}.spm.util.imcalc.input = cellstr(spm_select('expand',fullfile(sesd,'func',boldnii(funcs).name)));
                    matlabbatch{2}.spm.util.imcalc.output = [boldnii(funcs).name(1:end-4),'_var.nii'];
                    matlabbatch{2}.spm.util.imcalc.outdir = {fullfile(data_analysis_path,'QC',allsubs(sub).name)};
                    matlabbatch{2}.spm.util.imcalc.expression = 'var(X)';
                    matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
                    matlabbatch{2}.spm.util.imcalc.options.dmtx = 1;
                    matlabbatch{2}.spm.util.imcalc.options.mask = 0;
                    matlabbatch{2}.spm.util.imcalc.options.interp = 1;
                    matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
                    
                    spm_jobman('run', matlabbatch);
                    
                    delete([sesd,filesep,'func',filesep,boldniiunzip(funcs).name])
                end
            else
                fprintf('No func directory found.\n');
            end
        end
        
    end
end
fprintf('End\n')