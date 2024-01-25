% make script to fit diff models on start & end of trial 
% randomly sep in test & train segments
% test within chunks compared to between 

%% load in trial 
clear
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

cd('C:\Code\EphysCode'); 

if iscell(processed_behaviourData)
    processed_behaviourData = processed_behaviourData{1};
end

if iscell(processed_trialData)
    processed_trialData = processed_trialData{1};
end

%% break up in CL & OL segements 

 startOL = 178.189*1000; % sec
    endOL = 203.249*1000;
    lengthCL = 179.99*1000;
    
[CL, OL, CL_startStopIdx, OL_startStopIdx] = separateCL_OL(startOL, endOL, lengthCL, processed_behaviourData, processed_trialData);
%% break up into chunks w/ overlap

CL_chunks = [];
window = 60; 
window = window * 1000;
step = 30 * 1000;
count = 1; 
for chunk = 1:length(CL_startStopIdx)
    chunk_idx = CL_startStopIdx(chunk, 1):CL_startStopIdx(chunk, 2);
    for i = 1:step:length(chunk_idx) - window + 1
        CL_chunks(count,1) = chunk_idx(i);
        CL_chunks(count,2) = chunk_idx(i + window - 1); 
        count = count + 1;
    end
end

%%

seg_noMeno = [1:5];
seg_Meno =  (length(CL_chunks) - 9) :(length(CL_chunks) - 5);

    
 [B_noMeno, aveDeg_noMeno, rho_noMeno, prefHead_noMeno, R2_noMeno] = linearReg_chunks(seg_noMeno, processed_behaviourData, processed_trialData, CL_chunks, fileName);
 [B_Meno, aveDeg_Meno, rho_Meno, prefHead_Meno, R2_Meno] = linearReg_chunks(seg_Meno, processed_behaviourData, processed_trialData, CL_chunks, fileName);



fig = figure(1);clf
chunks = 1:length(B_Meno);
name = strcat(date,', ',cell_num,', ',trial,' Linear Regression Summary Plot');
sgtitle(strcat(name));

g(1) = subplot(4,2,1);
    plot(chunks,B_noMeno(1,:),'-ok','MarkerFaceColor', 'k')
    hold on
    plot(chunks,B_Meno(1,:),'-or','MarkerFaceColor', 'r')
    ylabel('BTh Vm')
g(2) = subplot(4,2,2);
    plot(chunks,B_noMeno(2,:),'-ok','MarkerFaceColor','k')
    hold on
    plot(chunks,B_Meno(2,:),'-or','MarkerFaceColor','r')
    ylabel('BVf Vm')
g(3) = subplot(4,2,3);
    plot(chunks,B_noMeno(3,:),'-ok','MarkerFaceColor', 'k')
    hold on
    plot(chunks,B_Meno(3,:),'-or','MarkerFaceColor', 'r')
    ylabel('Bvy Vm')
g(4) = subplot(4,2,4);
    plot(chunks,aveDeg_noMeno,'-ok','MarkerFaceColor', 'k')
    hold on
    plot(chunks,aveDeg_Meno,'-or','MarkerFaceColor', 'r')
    ylabel('heading')
    ylim([-200 200])
g(5) = subplot(4,2,5);
    plot(chunks,rho_noMeno,'-ok','MarkerFaceColor', 'k')
    hold on
    plot(chunks,rho_Meno,'-or','MarkerFaceColor', 'r')
    ylabel('vector str')
    ylim([0 1])
g(6) = subplot(4,2,6);
    plot(chunks,prefHead_noMeno,'-ok','MarkerFaceColor','k')
    hold on
    plot(chunks,prefHead_Meno,'-or','MarkerFaceColor','r')
    ylim([-200 200])
    ylabel('pHead Vm')  
g(7) = subplot(4,2,7); 
    plot(chunks,R2_noMeno,'-ok','MarkerFaceColor','k')
    hold on
    plot(chunks,R2_Meno,'-or','MarkerFaceColor','r')
    ylabel('R^2 Vm')
    text(chunks(end)-1,max(R2_Vm),num2str(mean(R2_Vm)))
    ylim([0 0.6])

linkaxes(g,'x')

