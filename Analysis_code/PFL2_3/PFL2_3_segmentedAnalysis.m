%% look at PFL2/3 tuning over specific parts of the trial
%clear
% rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys'; % ;'/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys' 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'%change dep on comp
% date = input('Date? ','s');
% cell_num = input('Cell? ','s');
% cell_num = strcat('cell_',cell_num);
% trial = input('Trial? ','s');
% trial = strcat('trial_',trial);
% fileName = fullfile(rootPath,date,cell_num,trial);

fileName = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL3\CL_OL\070121_1_R'; 
cd(fileName)
load('pro_trialData.mat');
load('pro_behaviourData.mat');
load('trialMeta.mat');
load('trialData.mat');

if isfield(trialMeta, 'notes')
    disp(trialMeta.notes)
end

cd('C:\Code\EphysCode');

%% remove parts of the trial where the fly is not moving
% must reload data before calling this section again 

% Change if only part of the trial is usable
start = 0;
finish = processed_behaviourData.time(end);

iStart = find(processed_behaviourData.time == start); 
iEnd = find(processed_behaviourData.time == finish); 

[no0Vel_frag, no0Vel_idx] = remove0velocity(1.5, 100, processed_behaviourData);

% simply gather all timepoints during which the fly's speed was above the
% threshold

       f = fieldnames(processed_behaviourData);
    for k=1:numel(f)
        no0vel_CL.(f{k}) = processed_behaviourData.(f{k})(no0Vel_idx); 
    end
    
           f = fieldnames(processed_trialData);
    for k=1:numel(f)
        no0vel_CL.(f{k}) = processed_trialData.(f{k})(no0Vel_idx); 
    end




%% recompute heatmaps to find pref angle

tStart = 1;
tEnd = length(behaviourData); 

iStart = 1; 
iEnd = length(behaviourData{1, 1}.vel_for); 
% %%
% Activity_Heatmap(tStart, tEnd, behaviourData, trialData, fileName, iStart, iEnd)
% 
% %% load in old heatmaps & line plots that included 0 vel segments and compare
% cd(fileName)
% openfig('directionPref_vel.fig');
% openfig('ActivityVsVelocity_splitHead.fig')
% cd('C:\Code\EphysCode');

%% recompute simple velocity/activity line plots

prefHead =-35; 
range = 120;
step = 0.5;
tStart = 1;
tEnd = length(processed_behaviourData); 

PFL2_3_behaviourVSactivity_lineplots(tStart, tEnd, prefHead, step, range, date, cell_num, trial, trialData{1}, behaviourData{1}, fileName, iStart, iEnd) 

%% recompute lagged line plots

% lag = -435; % ms
% t = 1; 
% prefHead = -80; 
% range = 100;
% step = 0.5;
% 
% PFL2_3_behaviourVSactivity_lag(t, lag, prefHead, step, range, date, cell_num, trial, trialData, behaviourData, fileName, iStart, iEnd) 

%% calculate average cross correlation across larger fragments

%% Segment trials into periods with different estimated menotaxic headings