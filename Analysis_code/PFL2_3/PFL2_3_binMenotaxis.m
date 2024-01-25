% load in exp data
close all
clear
rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys' ; % ;  '/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys'%change dep on comp
date = input('Date? ','s');
cell_num = input('Cell? ','s');
cell_num = strcat('cell_',cell_num);
trial = input('Trial? ','s');
trial = strcat('trial_',trial);
fileName = fullfile(rootPath,date,cell_num,trial);

cd(fileName)
load('pro_trialData.mat');
load('pro_behaviourData.mat');
openfig('directionPref_vel.fig');


cd('C:\Code\EphysCode\Analysis_code'); %'/Users/elenawesteinde/Documents/EphysCode'

%% Extract mean heading vectors & the windows in which they were calculated
window = 60;
minVel = 1.5; 
sampRate = 1000; 

[rho, theta, idx_windows] = PFL2_3_headingDist(window, minVel, processed_behaviourData, sampRate,fileName);
% rho = vector str
% theta = heading angle

%% bin data by their vector strength

% this is likely better done with some fancier clustering technique but for a
% quick & dirty test just use a set threshold to define menotaxis

% set the + & - menotaxis threshold, above + = menotaxis, below - = not &
% extract continuous segements where the fly is likely menotaxing

pos_meno = 0.6;
neg_meno = 0.4;

[pos_meno_frag, pos_meno_idx] = ThresholdData(pos_meno, 0, rho, 1, processed_trialData);
[neg_meno_Thres_frag, neg_meno_idx] = ThresholdData(neg_meno, 1, rho, 1, processed_trialData);

%% Collect all timepoints when fly was menotaxing & when it wasn't
pos_idx = [];
for i = 1:length(pos_meno_idx)
    start = (pos_meno_idx(i)-1)*window*sampRate; 
    finish = (pos_meno_idx(i)*window*sampRate)-1;
    pos_idx = [pos_idx, [start:finish]]; 
end
pos_idx = pos_idx + 1; %matlab indexing starts at 1

neg_idx = [];
for i = 1:length(neg_meno_idx)
    start = (neg_meno_idx(i)-1)*window*sampRate; 
    finish = (neg_meno_idx(i)*window*sampRate)-1;
    neg_idx = [neg_idx, [start:finish]]; 
end
neg_idx = neg_idx + 1; %matlab indexing starts at 1

pos_meno_data.vel_for = processed_behaviourData{1, 1}.vel_for([pos_idx]);
pos_meno_data.vel_side = processed_behaviourData{1, 1}.vel_side([pos_idx]);
pos_meno_data.vel_yaw = processed_behaviourData{1, 1}.vel_yaw([pos_idx]);
pos_meno_data.angle = processed_behaviourData{1, 1}.angle([pos_idx]);

pos_meno_data.fRate_sec = processed_trialData{1, 1}.fRate_sec([pos_idx]);
pos_meno_data.smooth_Vm = processed_trialData{1, 1}.smooth_Vm([pos_idx]);

neg_meno_data.vel_for = processed_behaviourData{1, 1}.vel_for([neg_idx]);
neg_meno_data.vel_side = processed_behaviourData{1, 1}.vel_side([neg_idx]);
neg_meno_data.vel_yaw = processed_behaviourData{1, 1}.vel_yaw([neg_idx]);
neg_meno_data.angle = processed_behaviourData{1, 1}.angle([neg_idx]);

neg_meno_data.fRate_sec = processed_trialData{1, 1}.fRate_sec([neg_idx]);
neg_meno_data.smooth_Vm = processed_trialData{1, 1}.smooth_Vm([neg_idx]);

%% re run basic analyses over the pos & neg meno data & compare 

%% recompute heatmaps to find pref angle

tStart = 1;
tEnd = length(pos_meno_data); 

posiStart = 1; 
posiEnd = length(pos_meno_data.vel_for); 

negiStart = 1; 
negiEnd = length(neg_meno_data.vel_for); 
%%
Activity_Heatmap(tStart, tEnd, pos_meno_data, pos_meno_data, fileName, posiStart, posiEnd)
Activity_Heatmap(tStart, tEnd, neg_meno_data, neg_meno_data, fileName, negiStart, negiEnd)

%% load in old heatmaps & line plots that included 0 vel segments and compare
cd(fileName)
openfig('directionPref_vel.fig');
openfig('ActivityVsVelocity_splitHead.fig')
cd('C:\Code\EphysCode');

%% recompute simple velocity/activity line plots

prefHead =-35; 
range = 120;
step = 0.5;

PFL2_3_behaviourVSactivity_lineplots(tStart, tEnd, prefHead, step, range, date, cell_num, trial, pos_meno_data, pos_meno_data, fileName, posiStart, posiEnd) 
PFL2_3_behaviourVSactivity_lineplots(tStart, tEnd, prefHead, step, range, date, cell_num, trial, neg_meno_data, neg_meno_data, fileName, negiStart, negiEnd) 

%% recompute lagged line plots

lag = -435; % ms
t = 1; 
prefHead = -80; 
range = 100;
step = 0.5;

PFL2_3_behaviourVSactivity_lag(t, lag, prefHead, step, range, date, cell_num, trial, pos_meno_data, pos_meno_data, fileName, posiStart, posiEnd) 

%% recompute correlation plots comparing strong & weak menotaxis boughts

% first choose what minute(s) of data to feed to xcorr function
pos_idx = [];
for i = 14:16
    start = (i-1)*window*sampRate; 
    finish = (i*window*sampRate)-1;
    pos_idx = [pos_idx, [start:finish]]; 
end
pos_idx = pos_idx + 1; %matlab indexing starts at 1

neg_idx = [];
for i = 18:20
    start = (i-1)*window*sampRate; 
    finish = (i*window*sampRate)-1;
    neg_idx = [neg_idx, [start:finish]]; 
end
neg_idx = neg_idx + 1; %matlab indexing starts at 1

pos_meno_corr.vel_for = processed_behaviourData{1, 1}.vel_for([pos_idx]);
pos_meno_corr.vel_side = processed_behaviourData{1, 1}.vel_side([pos_idx]);
pos_meno_corr.vel_yaw = processed_behaviourData{1, 1}.vel_yaw([pos_idx]);
pos_meno_corr.angle = processed_behaviourData{1, 1}.angle([pos_idx]);

pos_meno_corr.fRate_sec = processed_trialData{1, 1}.fRate_sec([pos_idx]);
pos_meno_corr.smooth_Vm = processed_trialData{1, 1}.smooth_Vm([pos_idx]);

neg_meno_corr.vel_for = processed_behaviourData{1, 1}.vel_for([neg_idx]);
neg_meno_corr.vel_side = processed_behaviourData{1, 1}.vel_side([neg_idx]);
neg_meno_corr.vel_yaw = processed_behaviourData{1, 1}.vel_yaw([neg_idx]);
neg_meno_corr.angle = processed_behaviourData{1, 1}.angle([neg_idx]);

neg_meno_corr.fRate_sec = processed_trialData{1, 1}.fRate_sec([neg_idx]);
neg_meno_corr.smooth_Vm = processed_trialData{1, 1}.smooth_Vm([neg_idx]);


prefHead = -80; 
range = 120;

posiStart = 1; 
posiEnd = length(pos_meno_corr.vel_for); 

negiStart = 1; 
negiEnd = length(neg_meno_corr.vel_for); 

xcorr_headingBin(prefHead,range, pos_meno_corr, pos_meno_corr, fileName,posiStart,posiEnd)
xcorr_headingBin(prefHead,range, neg_meno_corr, neg_meno_corr, fileName,negiStart,negiEnd)


%% remove parts of the trial where the fly's speed is below some threshold
% must reload data before calling this section again 

% Change if only part of the trial is usable

[pos_meno_no0Vel_frag, pos_meno_no0Vel_idx] = remove0velocity(1.5, 100, pos_meno_data, pos_meno_data);
[neg_meno_no0Vel_frag, neg_meno_no0Vel_idx] = remove0velocity(1.5, 100, neg_meno_data, neg_meno_data);

% simply gather all timepoints during which the fly's speed was above the
% threshold

pos_meno_no0Vel.vel_for = pos_meno_data.vel_for(pos_meno_no0Vel_idx);
pos_meno_no0Vel.vel_yaw = pos_meno_data.vel_yaw(pos_meno_no0Vel_idx);
pos_meno_no0Vel.vel_side = pos_meno_data.vel_side(pos_meno_no0Vel_idx);
pos_meno_no0Vel.angle = pos_meno_data.angle(pos_meno_no0Vel_idx);
pos_meno_no0Vel.smooth_Vm = pos_meno_data.smooth_Vm(pos_meno_no0Vel_idx);
pos_meno_no0Vel.fRate_sec = pos_meno_data.fRate_sec(pos_meno_no0Vel_idx);

neg_meno_no0Vel.vel_for = neg_meno_data.vel_for(neg_meno_no0Vel_idx);
neg_meno_no0Vel.vel_yaw = neg_meno_data.vel_yaw(neg_meno_no0Vel_idx);
neg_meno_no0Vel.vel_side = neg_meno_data.vel_side(neg_meno_no0Vel_idx);
neg_meno_no0Vel.angle = neg_meno_data.angle(neg_meno_no0Vel_idx);
neg_meno_no0Vel.smooth_Vm = neg_meno_data.smooth_Vm(neg_meno_no0Vel_idx);
neg_meno_no0Vel.fRate_sec = neg_meno_data.fRate_sec(neg_meno_no0Vel_idx);

posiStart_no0vel = 1; 
posiEnd_no0vel = length(pos_meno_no0Vel.vel_for); 

negiStart_no0vel = 1; 
negiEnd_no0vel = length(neg_meno_no0Vel.vel_for); 

%% recompute line plots with timepoints of 0 speed removed

prefHead =-80; 
range = 120;
step = 0.5;
tStart = 1;
tEnd = length(pos_meno_data); 

PFL2_3_behaviourVSactivity_lineplots(tStart, tEnd, prefHead, step, range, date, cell_num, trial, pos_meno_no0Vel, pos_meno_no0Vel, fileName, posiStart_no0vel, posiEnd_no0vel) 
PFL2_3_behaviourVSactivity_lineplots(tStart, tEnd, prefHead, step, range, date, cell_num, trial, neg_meno_no0Vel, neg_meno_no0Vel, fileName, negiStart_no0vel, negiEnd_no0vel) 









