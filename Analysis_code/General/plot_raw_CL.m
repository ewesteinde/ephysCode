rootPath = 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'; %'/Users/ewesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys' ; %change dep on comp
date = input('Date? ','s');
cell_num = input('Cell? ','s');
cell_num = strcat('cell_',cell_num);
trial = input('Trial? ','s');
trial = strcat('trial_',trial);
fileName = fullfile(rootPath,date,cell_num,trial);

cd(fileName)

%load('rawData.mat');
load('trialData.mat');
load('trialMeta.mat');
load('behaviorData.mat'); 

%%
figure(t); clf;
    h(1) = subplot(8,1,1);
    plot(trialData{t}.time, trialData{t}.current, 'k')
    ylabel('Current (pA)')

    h(2) = subplot(8,1,2:4);
    plot(trialData{t}.time, trialData{t}.scaledOutput, 'k')
    ylabel('Voltage (mV)')
    xlabel('Time (s)')

    sgtitle(['Trial ' num2str(t)])
    hold on 
           
    h(3) = subplot(8,1,5);
    plot(behaviorData{t}.time, behaviorData{t}.angle, 'k')
    ylabel('Pattern Angle')

    h(4) = subplot(8,1,6);
    plot(behaviorData{t}.time, behaviorData{t}.vel_for, 'k')
    ylabel('Vel For')
    
    h(5) = subplot(8,1,7);
    plot(behaviorData{t}.time, behaviorData{t}.vel_yaw, 'k')
    ylabel('Vel Yaw')
    
    h(6) = subplot(8,1,8);
    plot(behaviorData{t}.time, behaviorData{t}.vel_side, 'k')
    ylabel('Vel Side')
    xlabel('Time (s)')

    linkaxes(h,'x');
    
 