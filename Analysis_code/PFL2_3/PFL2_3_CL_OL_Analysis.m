clear
rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL3\CL_OL'; % ;'/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys' 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'%change dep on comp
 % ;'/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys' 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'%change dep on comp

date = input('Date? ','s');
cell_num = input('Cell? ','s');
cell_num = strcat('cell_',cell_num);
trial = input('Trial? ','s');
trial = strcat('trial_',trial);
%fileName = fullfile(rootPath,date,cell_num,trial);

cd(fileName)
load('pro_trialData.mat');
load('pro_behaviourData.mat');
load('trialMeta.mat');
load('trialData.mat');

if isfield(trialMeta, 'notes')
    disp(trialMeta.notes)
end

cd('C:\Code\EphysCode'); 

if iscell(processed_trialData)
    processed_trialData = processed_trialData{1}; 
end

if iscell(processed_behaviourData)
     processed_behaviourData = processed_behaviourData{1}; 
end


%% extract signals for signal analyzer

% vf = processed_behaviourData.vel_for;
% vf(isnan(vf)) = 0;
% vs = processed_behaviourData.vel_side;
% vs(isnan(vs)) = 0;
% processed_behaviourData.vel_side(isnan(processed_behaviourData.vel_side)) = 0;
% vy = processed_behaviourData.vel_yaw;
% vy(isnan(vy)) = 0;
% processed_behaviourData.vel_yaw(isnan(processed_behaviourData.vel_yaw)) = 0;
% angle = processed_behaviourData.angle; 
% speed = sqrt(vf.^2 + vs.^2); 
% 
% ScO = processed_trialData.scaledOutput_down; 
% 
% Vm = processed_trialData.smooth_Vm; 
% 
% %%
plotAngle_overlayVel_ScO(processed_trialData, processed_behaviourData)
        
%% Change if only part of the trial is usable
start =0;
finish = processed_behaviourData.time(end);

iStart = find(processed_behaviourData.time == start); 
iEnd = find(processed_behaviourData.time == finish); 


f = fieldnames(processed_behaviourData);
for k=1:numel(f)
    processed_behaviourData.(f{k}) = processed_behaviourData.(f{k})(iStart:iEnd); 
end


f = fieldnames(processed_trialData);
for k=1:numel(f)
    processed_trialData.(f{k}) = processed_trialData.(f{k})(iStart:iEnd); 
end
        
    %% separate idx at which fly is moving in CL vs OL
    
    startOL = 179.442*1000; % sec
    endOL = 204.521*1000;
    lengthCL = 179.99*1000;
    
% start = 1819.45;
% finish = 1844.47;
%     startOL  = find(processed_behaviourData.time == start); 
%     endOL = find(processed_behaviourData.time == finish); 

    [CL, OL, CL_startStopIdx] = separateCL_OL(startOL, endOL, lengthCL, processed_behaviourData, processed_trialData);

    save(fullfile(fileName,'CL_idx.mat'),'CL_startStopIdx')
    save(fullfile(fileName,'CL.mat'),'CL')
    
    
%% run standard analyses on CL & OL segments

%% break into diff chunks if wanted


% heading preference heatmaps & line plots
tStart = 1;
tEnd = 1;

prefHead = 0; 

Activity_Heatmap(10, tStart, tEnd, CL, CL, fileName)
Activity_Heatmap(30, tStart, tEnd, OL, OL, fileName)


figure();clf;
Vm = CL.smooth_Vm;
angle = CL.angle;
edges = [-180:30:180];
[centers, mean_bin] = create_binned_mean(angle, Vm, edges);
plot(centers, mean_bin,'-o');
xlabel('Angle');
ylabel('Vm'); 
title('CL heading pref')

keep = input('Save? ','s');

if strcmp(keep, 'y')
    cd(fileName) 
    saveas(gcf,'CL_HeadingPref_lineplot.fig')
    cd('C:\Code\EphysCode\Analysis_code') %'/Users/elenawesteinde/Documents/EphysCode/Analysis_code' 'C:\Code\EphysCode'
end



figure(4);clf;
Vm = OL.smooth_Vm;
angle = OL.angle;
edges = [-180:30:180];
[centers, mean_bin] = create_binned_mean(angle, Vm, edges);
plot(centers, mean_bin,'-o');
xlabel('Angle');
ylabel('Vm'); 
title('OL heading pref')

keep = input('Save? ','s');

if strcmp(keep, 'y')
    cd(fileName) 
    saveas(gcf,'OL_HeadingPref_lineplot.fig')
    cd('C:\Code\EphysCode\Analysis_code') %'/Users/elenawesteinde/Documents/EphysCode/Analysis_code' 'C:\Code\EphysCode'
end


%% activity-velocity line plots

CL_all = CL; 
CL_goalIdx = 1:254676;%:length(CL_all.vel_for);

fn = fieldnames(CL_all);
for k=1:numel(fn)
        CL.(fn{k}) = CL_all.(fn{k})(CL_goalIdx);
end

prefHead = -20; 
range = 120;
step = 20;
tStart = 1;
tEnd = 1; 
PFL2_3_behaviourVSactivity_lineplots_poster(prefHead, step, range, CL, CL) ;


PFL2_3_behaviourVSactivity_lineplots(prefHead, step, range, date, cell_num, trial, CL, CL) ;


keep = input('Save? ','s');

if strcmp(keep, 'y')
    cd(fileName) 
    savefig(a,'ActivityVsVelocity_splitHead_CL.fig')
    savefig(b,'ActivityVsVelocity_allHead_CL.fig')
    cd('C:\Code\EphysCode\Analysis_code') %'C:\Code\EphysCode'
end
[c, d] = PFL2_3_behaviourVSactivity_lineplots(tStart, tEnd, prefHead, step, range, date, cell_num, trial, OL, OL, fileName) ;

keep = input('Save? ','s');

if strcmp(keep, 'y')
    cd(fileName) 
    savefig(c,'ActivityVsVelocity_splitHead_OL.fig')
    savefig(d,'ActivityVsVelocity_allHead_OL.fig')
    cd('C:\Code\EphysCode\Analysis_code') %'C:\Code\EphysCode'
end

repeatno0vel = input('repeat line plots with timepoints of no movement removed? y/n ','s');
if strcmp(repeatno0vel,'y')
    [no0Vel_frag, no0Vel_idx] = remove0velocity(1.5, 0,  CL);

    % simply gather all timepoints during which the fly's speed was above the
    % threshold

       f = fieldnames(CL);
    for k=1:numel(f)
        no0vel_CL.(f{k}) = CL.(f{k})(no0Vel_idx); 
    end

    Activity_Heatmap(10, tStart, tEnd, no0vel_CL, no0vel_CL, fileName)
    
    [e,f] = PFL2_3_behaviourVSactivity_lineplots(tStart, tEnd, prefHead, step, range, date, cell_num, trial, no0vel_CL, no0vel_CL, fileName);

    keep = input('Save? ','s');

    if strcmp(keep, 'y')
        cd(fileName) 
        savefig(e,'ActivityVsVelocity_splitHead_CL_no0vel.fig') 
        savefig(f,'ActivityVsVelocity_allHead_CL_no0vel.fig')
        cd('C:\Code\EphysCode\Analysis_code') %'C:\Code\EphysCode'
    end
    
    window = 60; %seconds
    minVel = 2; %mm/s
    sampRate = 1000; %Hz 


    [rho, theta, idx_windows, g, h] = PFL2_3_headingDist(window, minVel, no0vel_CL, sampRate, fileName, 1);
    
    keep = input('Save? ','s');

    if strcmp(keep, 'y')
        cd(fileName) 
        savefig(g,strcat('no0vel_HeadingDist ',num2str(window/1000),'s window.fig'))
        savefig(h,strcat('no0vel_HeadingDistHist ',num2str(window/1000),'s window.fig')) 
        cd('C:\Code\EphysCode') %'C:\Code\EphysCode'
    end
end

%% plot heading distributions during CL segements of the fly across the trial 
window = 60; %seconds
minVel = 2; %mm/s
sampRate = 1000; %Hz 


[rho, theta, idx_windows, i, k] = PFL2_3_headingDist(window, minVel, CL, sampRate, fileName,1);

keep = input('Save? ','s');

    if strcmp(keep, 'y')
        cd(fileName) 
        savefig(i,strcat('HeadingDist ',num2str(window/1000),'s window.fig'))
        savefig(k,strcat('HeadingDistHist ',num2str(window/1000),'s window.fig')) 
        cd('C:\Code\EphysCode') %'C:\Code\EphysCode'
    end

%% Calculate cross correlation

basic_xcorr_PFL( processed_trialData, processed_behaviourData, date, cell_num, trial, fileName);

basic_xcorr_PFL( CL, CL, date, cell_num, trial, fileName);


%% plot activity-behaviour line plots across headings & at different lags
% negative lag means neural activity precedes behaviour, positive means
% neural activity follows behaviour 

lag = 150; % ms
t = 1; 
prefHead = 65; 
range = 100;
step = 0.5;

PFL2_3_behaviourVSactivity_lag(t, lag, prefHead, step, range, date, cell_num, trial, CL, CL, fileName, iStart_CL, iEnd_CL) 
PFL2_3_behaviourVSactivity_lag(t, lag, prefHead, step, range, date, cell_num, trial, OL, OL, fileName, iStart_OL, iEnd_OL) 




