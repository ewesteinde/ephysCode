clear
close all 

rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys'; % ;'/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys' 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'%change dep on comp
 % ;'/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys' 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'%change dep on comp

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

cd('C:\Code\EphysCode'); 

%% extract signals for signal analyzer

vf = processed_behaviourData.vel_for;
vf(isnan(vf)) = 0;
vs = processed_behaviourData.vel_side;
vs(isnan(vs)) = 0;
processed_behaviourData.vel_side(isnan(processed_behaviourData.vel_side)) = 0;
vy = processed_behaviourData.vel_yaw;
vy(isnan(vy)) = 0;
processed_behaviourData.vel_yaw(isnan(processed_behaviourData.vel_yaw)) = 0;
angle = processed_behaviourData.angle; 
speed = sqrt(vf.^2 + vs.^2); 

ScO = processed_trialData.scaledOutput_down; 

Vm = processed_trialData.smooth_Vm; 



%% visualize individual trials

t = 1;

%
nActivity = processed_trialData{t}.scaledOutput_down;
figure();clf;      
        
        h(1) = subplot(4,1,1);
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.angle, 'k') 
        ylabel('angle')
       
        
        h(2) = subplot(4,1,2);
        yyaxis left
        plot(processed_behaviourData{t}.time,nActivity, 'k') 
        ylabel('Vm')
        %ylim([-70 -45])
        hold on 
        yyaxis right
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_for, 'r')
        ylim([-(max(processed_behaviourData{t}.vel_for)) max(processed_behaviourData{t}.vel_for)])
        ylabel('Vf mm/sec')
        
        
        
        h(3) = subplot(4,1,3);
        yyaxis left
        plot(processed_behaviourData{t}.time,nActivity, 'k')
        ylabel('Vm')
       % ylim([-70 -45])
        hold on 
        yyaxis right
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_side, 'r')
        ylim([-(max(processed_behaviourData{t}.vel_side)) max(processed_behaviourData{t}.vel_side)])
        ylabel('Vs mm/sec')
        
        h(4) = subplot(4,1,4);
        yyaxis left
        plot(processed_behaviourData{t}.time,nActivity, 'k')
        ylabel('Vm')

        %ylim([-55 -])
        hold on 
        yyaxis right
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_yaw, 'r')
        ylim([(min(processed_behaviourData{t}.vel_yaw)) max(processed_behaviourData{t}.vel_yaw)])
        ylabel('Vy deg/sec')
        
        
        linkaxes(h,'x');


%% Change if only part of the trial is usable
start = 0;
finish = processed_behaviourData.time(end);

iStart = find(processed_behaviourData.time == start); 
iEnd = find(processed_behaviourData.time == finish); 

%% use to remove slow linear depolarization in Vm over the trial period

dt_smooth_Vm = detrend(processed_trialData.smooth_Vm);
dt_Vm = detrend(processed_trialData.Vm);
baseValue = median(processed_trialData.smooth_Vm(1:(10*100)));

figure()
subplot(3,1,1)
plot(processed_trialData{1}.smooth_Vm)
subplot(3,1,2)
plot(dt_smooth_Vm)
subplot(3,1,3)
plot(dt_smooth_Vm + baseValue)

processed_trialData{1}.smooth_Vm = dt_smooth_Vm + baseValue;
processed_trialData{1}.Vm = dt_Vm;

%% Plot heatmaps b/w neural activity & velocity/speed
%close all 

tStart = 1;
tEnd = 1; 
angleBin = 10;  
Activity_Heatmap(angleBin, tStart, tEnd, processed_behaviourData, processed_trialData, fileName)

%% plot line plots across headings & compared from estimated pref vs not pref headings

cd(fileName)
openfig('directionPref_vel.fig');
cd('/Users/elenawesteinde/Documents/EphysCode/Analysis_code');
%%

prefHead = 150; 
range = 120;
step = 0.5;
tStart = 1;
tEnd = length(processed_behaviourData); 

PFL2_3_behaviourVSactivity_lineplots(tStart, tEnd, prefHead, step, range, date, cell_num, trial, processed_trialData, processed_behaviourData, fileName, iStart, iEnd) 

%% Calculate cross correlation

basic_xcorr_PFL(processed_trialData, processed_behaviourData, date, cell_num, trial, fileName)

%% plot activity-behaviour line plots across headings & at different lags
% negative lag means neural activity precedes behaviour, positive means
% neural activity follows behaviour 


lag = 100; % ms
t = 1; 
prefHead = 65; 
range = 100;
step = 0.5;

PFL2_3_behaviourVSactivity_lag(t, lag, prefHead, step, range, date, cell_num, trial, processed_trialData, processed_behaviourData, fileName, iStart, iEnd) 

%% plot heading distributions of the fly across the trial 

window = 60; %seconds
minVel = 2; %mm/s
sampRate = 1000; %Hz  

PFL2_3_headingDist(window, minVel, processed_behaviourData, sampRate, fileName);

%% calculate cross correlation binned by heading

start = 0;
finish = processed_behaviourData{1}.time(end);

prefHead = 65; 
range = 100;

iStart = find(processed_behaviourData{1}.time == start); 
iEnd = find(processed_behaviourData{1}.time == finish); 

xcorr_headingBin(prefHead,range, processed_trialData{1}, processed_behaviourData{1}, fileName,iStart,iEnd)








 



