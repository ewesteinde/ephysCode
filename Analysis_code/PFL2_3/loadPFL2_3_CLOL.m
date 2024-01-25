function [processed_trialData, processed_behaviourData, fileName] = loadPFL2_3_CLOL(cellType, folderName)
% loads in CL_OL PFL2/3 data

if ~exist('folderName', 'var')
    if strcmp(cellType,'PFL3')
        rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL3\CL_OL'; 
        date = input('Date? ','s');
        cell_num = input('Cell? ','s');
        LAL = input('LAL L/R? ','s');
        fileName = fullfile(rootPath,strcat(date,'_',cell_num,'_',LAL));
    elseif strcmp(cellType,'PFL2')
        rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL2\CL_OL'; 
        date = input('Date? ','s');
        cell_num = input('Cell? ','s');
        fileName = fullfile(rootPath,strcat(date,'_',cell_num));
    end
else 
    if strcmp(cellType,'PFL3')
        rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL3\CL_OL'; 
    elseif strcmp(cellType,'PFL2')
        rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL2\CL_OL'; 
    end
    fileName = fullfile(rootPath,folderName);
end

cd(fileName)
load('pro_trialData.mat');
load('pro_behaviourData.mat');
load('trialMeta.mat');

if iscell(processed_behaviourData)
    processed_behaviourData = processed_behaviourData{1};
end

if iscell(processed_trialData)
    processed_trialData = processed_trialData{1};
end

cd('C:\Code\EphysCode');