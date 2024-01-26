function [processed_trialData, processed_behaviourData, fileName] = loadPFL2_3_CLJump(cellType, folderName)
% loads in earlier PFL2/3 data prior to when I started doing CL_OL
% experiments & cocatinates cells in trials where I had 3 consecutive
% trials split into 3 cells. Don't use for cross correlation. 
if ~exist('folderName', 'var')
    if strcmp(cellType,'PFL3')
        rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL3\CL_jump'; 
        date = input('Date? ','s');
        cell_num = input('Cell? ','s');
        LAL = input('LAL L/R? ','s');
        fileName = fullfile(rootPath,strcat(date,'_',cell_num,'_',LAL));
    elseif strcmp(cellType,'PFL2')
        rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL2\CL_jump'; 
        date = input('Date? ','s');
        cell_num = input('Cell? ','s');
        fileName = fullfile(rootPath,strcat(date,'_',cell_num));
    end
else 
    if strcmp(cellType,'PFL3')
        rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL3\CL_jump'; 
    elseif strcmp(cellType,'PFL2')
        rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL2\CL_jump'; 
    end
    fileName = fullfile(rootPath,folderName);
end

cd(fileName)
load('pro_trialData.mat');
load('pro_behaviourData.mat');
load('trialMeta.mat');

cd('C:\Code\EphysCode'); 

if iscell(processed_behaviourData) && length((processed_behaviourData)) > 1 
    processed_behaviourData_temp = [];
    processed_trialData_temp = [];
    
    f = fieldnames(processed_behaviourData{1});
    for k=1:numel(f)
        processed_behaviourData_temp.(f{k}) = []; 
    end

    f = fieldnames(processed_trialData{1});
    for k=1:numel(f)
        processed_trialData_temp.(f{k}) = []; 
    end
    
    for c = 1:length(processed_behaviourData)
        f = fieldnames(processed_behaviourData{1});
        for k=1:numel(f)
            processed_behaviourData_temp.(f{k}) = [processed_behaviourData_temp.(f{k}) processed_behaviourData{c}.(f{k})]; 
        end
        
        f = fieldnames(processed_trialData{1});
        for k=1:numel(f)
            processed_trialData_temp.(f{k}) = [processed_trialData_temp.(f{k}) processed_trialData{c}.(f{k})]; 
        end
    end
    
    processed_trialData = processed_trialData_temp;
    processed_behaviourData = processed_behaviourData_temp;
  
elseif iscell(processed_behaviourData) && iscell(processed_trialData)
    processed_behaviourData = processed_behaviourData{1};
    processed_trialData = processed_trialData{1};
    
end