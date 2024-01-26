%% load in experiment 

clear
rootPath = '/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys'; % ;'/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys' 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'%change dep on comp
date = input('Date? ','s');
cell_num = input('Cell? ','s');
cell_num = strcat('cell_',cell_num);
trial = input('Trial? ','s');
trial = strcat('trial_',trial);
fileName = fullfile(rootPath,date,cell_num,trial);

cd(fileName)
load('pro_trialData.mat');
load('pro_behaviourData.mat');
load('trialMeta.mat');
load('trialData.mat');

if isfield(trialMeta, 'notes')
    disp(trialMeta.notes)
end

cd('/Users/elenawesteinde/Documents/EphysCode'); %  'C:\Code\EphysCode' '/Users/elenawesteinde/Documents/EphysCode'%change dep on comp



%% compute translation velocity in preferred direction of the cell
offset = 31;
[processed_behaviourData] = calc_velDot(offset, processed_trialData, processed_behaviourData, fileName);

%% Calculate percentage of trial during which the fly was moving
close all 

minVel = 0.5;
total_above_vel = 0; 
total_trial_length = 0; 

for t = 1:length(processed_trialData)
    speed_singleTrial = abs(processed_behaviourData{t}.vel_dot); 
    above_vel_min = speed_singleTrial(speed_singleTrial > minVel); 
    num_aboveVelMin = length(above_vel_min); 
    total_above_vel = total_above_vel + num_aboveVelMin; 
    total_trial_length = total_trial_length + length(processed_behaviourData{t}.vel_dot); 
end 

percent_above_vel = (total_above_vel/total_trial_length) * 100;



%% calculate cross correlation b/w neural activity & translation velocity
close all 
[corrData] = basic_xcorr(trialData, processed_trialData, processed_behaviourData, date, cell_num, trial, fileName); 

%% estimate preferred heading direction

VelDot_Heatmap(1, 3, processed_behaviourData, processed_trialData, fileName);

%% plot neural activity - velocity/speed relationships dependent on heading
close all 
cd(fileName)
openfig('ActivityVsVelocity.fig');
openfig('directionPref_velDot.fig');
cd('/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys');
%%
close all 

Tstart = 1;
Tend = 3; 
prefHead = 60;
offset = -31;

step = 0.5; 
range = 120; 
maxValPlot = 6;
minValPlot = -2;


[saveBinsFRVel, saveBinsVmVel] = plotActivityvsSpeed_prefDirection(Tstart, Tend, prefHead, step, range, offset, maxValPlot, minValPlot, date, cell_num, trial, processed_trialData, processed_behaviourData,fileName) ;
%%
trendSum{count}{1} = saveBinsVmVel;
trendSum{count}{2} = saveBinsFRVel;
count = count + 1; 







