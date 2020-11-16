%% Cleanup IntendedFor entries in fieldmap-accompanying json-sidecar file
% json-files from fieldmaps with IntendedFor entries for tasks which do not
% exist in a session prevent BIDS validation. RaBIDS v0.2 writes tasks to IntendedFor fields, regardless whether they exist or not.
% Use this program to remove tasks that do not exist in a session.

% Christian Paret, ZI Mannheim, 2020

clear

%% Select dataset
% Use below routine if dataset has been created with RaBIDS and run it from the code directory. Otherwise: define parent directory of your BIDS dataset in variable 'data_analysis_path', add path to saveJSONfile, and comment out lines 12-21.
data = readtable('datasheet.xlsx','ReadRowNames',true,'PreserveVariableNames',true,'NumHeaderLines',0);
AddPathline = find(contains(data.Properties.RowNames,'add path')); % add json path
userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));
for i = 1:length(AddPathline)
    addpath(data{AddPathline(i),userInputcol}{:});
end
userInputcol = find(strcmp(data.Properties.VariableNames,'UserInput'));
DataAnalysisPathline = find(strcmp(data.Properties.RowNames,'data analysis path'));
data_analysis_path = data{DataAnalysisPathline,userInputcol}{:};
datasetd = [data_analysis_path,filesep,'dataset'];
%%

allsubs = dir([datasetd,filesep,'sub-*']);

for sub = 1:length(allsubs)
    subd = [datasetd,filesep,allsubs(sub).name];
    
    if isfolder(subd)
        
        fprintf(['\n-------- ',allsubs(sub).name,' --------\n'])
        
        allses = dir([subd,filesep,'ses-*']);
        
        if ~isempty(allses)
            nses = length(allses);
            sesdir = 1;
        else
            allses(1).name = ['.',filesep];
            nses = 1;
            sesdir = 0;
        end
        
        for ses = 1:nses
            clear IFtaskid existtaskid
            sesd = [subd,filesep,allses(ses).name];
            if nses == 1
                fprintf('One scan session.\n')
            else
                fprintf([allses(ses).name,'\n'])
            end
            
            if isfolder([sesd,filesep,'fmap'])
                
                fmapjson = dir([sesd,filesep,'fmap',filesep,allsubs(sub).name,'_',allses(ses).name,'_phasediff.json']);
                
                if isempty(fmapjson)
                    fprintf('No json-sidecar found for fieldmap phasedifference image.\n');
                else
                    
                    % decode fmap json file
                    jsonfname = [sesd,filesep,'fmap',filesep,allsubs(sub).name,'_',allses(ses).name,'_phasediff.json'];
                    jsonf = jsondecode(fileread(jsonfname));
                    
                    if isfield(jsonf,{'IntendedFor'})
                        for bold = 1:length(jsonf.IntendedFor)
                            dum1 = strfind(jsonf.IntendedFor{bold},'task-');
                            dum2 = strfind(jsonf.IntendedFor{bold},'_bold');
                            IFtaskid{bold} = jsonf.IntendedFor{bold}(dum1+5:dum2-1);
                        end
                        
                        jsonf = rmfield(jsonf,{'IntendedFor'});
                        
                        % read existing tasks for this session from func directory
                        existtaskid = cell(1); % initiate variable
                        boldjson = dir([sesd,filesep,'func',filesep,allsubs(sub).name,'_',allses(ses).name,'_task-*_bold.json']);
                        if isempty(boldjson)
                            fprintf('No json-sidecar found for bold scans. Fieldmap not intended for any task.\n');
                        else
                            for bold = 1:length(boldjson)
                                dum1 = strfind(boldjson(bold).name,'task-');
                                dum2 = strfind(boldjson(bold).name,'_bold');
                                existtaskid{bold} = boldjson(bold).name(dum1+5:dum2-1);
                            end
                            
                            if length(IFtaskid)>length(existtaskid)
                                
                                % write existing tasks to fmap json
                                usetaskid = intersect(existtaskid,IFtaskid);
                                
                                if ~isempty(usetaskid)
                                    for bold = 1:length(usetaskid)
                                        if sesdir
                                            jsonf.IntendedFor{bold,1} = [allses(ses).name,'/func/',allsubs(sub).name,'_',allses(ses).name,'_task-',usetaskid{bold},'_bold.nii.gz'];
                                        else
                                            jsonf.IntendedFor{bold,1} = ['func/',allsubs(sub).name,'_',allses(ses).name,'_task-',usetaskid{bold},'_bold.nii.gz'];
                                        end
                                    end
                                end
                                fprintf('Remove non-existant tasks from field IntendedFor...\n');
                            else
                                fprintf('Bold scans exist for all tasks in field IntendedFor.\n');
                            end
                        end
                        saveJSONfile(jsonf,jsonfname) % Lior Kirsch (2020). Structure to JSON (https://www.mathworks.com/matlabcentral/fileexchange/50965-structure-to-json), MATLAB Central File Exchange. Retrieved February 27, 2020.
                    else
                        fprintf('Existing fieldmap not intended for any task.\n')
                    end
                end
                
            else
                fprintf('No fmap directory found.\n');
            end
        end
        
    end
end
fprintf('End\n')