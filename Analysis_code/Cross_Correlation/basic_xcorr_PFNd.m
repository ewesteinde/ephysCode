%% Calculate cross correlation b/w VelDot & Vm

% %load in desired trial
% %close all
% clear
% rootPath = '/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys'; % ; 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys' %change dep on comp
% date = input('Date? ','s');
% cell_num = input('Cell? ','s');
% cell_num = strcat('cell_',cell_num);
% trial = input('Trial? ','s');
% trial = strcat('trial_',trial);
% fileName = fullfile(rootPath,date,cell_num,trial);
% 
% cd(fileName)
% load('pro_trialData.mat');
% load('pro_behaviourData.mat');
% load('trialData.mat')
% load('trialMeta.mat')
% 
% cd('/Users/elenawesteinde/Documents/EphysCode'); %'C:\Code\EphysCode'

%% calculate & plot cross correlation for each trial
function [corrData] = basic_xcorr_PFNd(trialData, processed_trialData, processed_behaviourData, date, cell_num, trial, fileName) 

name = strcat(date,', ',cell_num,', ',trial);
      

ephysSettings;
downsample_Hz = 1000;  

    count = 1; 
    saveBinsVmVel = cell(1,3);
    saveBinsFRVel = cell(1,3); 
    saveBinsSOVel = cell(1,3);

for t = 1:length(processed_trialData)

%seemingly randomly some exps have a continuous stream of NaN values at
%start or end of trials, will need to cut these sections in order to run

velDot = processed_behaviourData{t}.vel_dot;  
smVm = processed_trialData{t}.smooth_Vm;
fRate = processed_trialData{t}.fRate_sec;
ScO = resample(trialData{t}.scaledOutput, downsample_Hz, settings.sampRate);
nanIDX  = find( isnan(velDot) );

if ~isempty(nanIDX)
    if nanIDX(1) == 1
        velDot = velDot(nanIDX(end)+1:length(velDot)); 
        smVm = smVm(nanIDX(end)+1:length(smVm)); 
        fRate = fRate(nanIDX(end)+1:length(fRate)); 
        ScO = ScO(nanIDX(end)+1:length(ScO)); 
    elseif nanIDX(end) == length(velDot)
        velDot = velDot(1:nanIDX(1)-1);
        smVm = smVm(1:nanIDX(1)-1); 
        fRate = fRate(1:nanIDX(1)-1); 
        ScO = ScO(1:nanIDX(1)-1); 
    else
        error('NaN values inside velDot')
    end
end
  
  
  
    [cVm, lagsVm] = xcorr(smVm,velDot,1000);
    [cFR, lagsFR] = xcorr(fRate,velDot,1000);
    [cSO, lagsSO] = xcorr(ScO,velDot,1000);

    saveBinsVmVel{count} = [(lagsVm./1000)', cVm];
    saveBinsFRVel{count} = [(lagsFR./1000)', cFR];
    saveBinsSOVel{count} = [(lagsSO./1000)', cSO];
    
    count = count + 1;
end

h = figure(20); clf; 

    subplot(3,3,1);
    plot(saveBinsSOVel{1}(:,1),saveBinsSOVel{1}(:,2))
    xline(saveBinsSOVel{1}(saveBinsSOVel{1}(:,2) == max(saveBinsSOVel{1}(:,2)),1),'-',num2str(saveBinsSOVel{1}(saveBinsSOVel{1}(:,2) == max(saveBinsSOVel{1}(:,2)),1)))
    legend('hide')
    title('T1 scaledOutput')
    subplot(3,3,4);
    plot(saveBinsSOVel{2}(:,1),saveBinsSOVel{2}(:,2))
    xline(saveBinsSOVel{2}(saveBinsSOVel{2}(:,2) == max(saveBinsSOVel{2}(:,2)),1),'-',num2str(saveBinsSOVel{2}(saveBinsSOVel{2}(:,2) == max(saveBinsSOVel{2}(:,2)),1)))
    legend('hide')
    title('T2 scaledOutput')
    subplot(3,3,7);
    plot(saveBinsSOVel{3}(:,1),saveBinsSOVel{3}(:,2))
    xline(saveBinsSOVel{3}(saveBinsSOVel{3}(:,2) == max(saveBinsSOVel{3}(:,2)),1),'-',num2str(saveBinsSOVel{3}(saveBinsSOVel{3}(:,2) == max(saveBinsSOVel{3}(:,2)),1)))
    legend('hide')
    title('T3 scaledOutput')

    subplot(3,3,2);
    plot(saveBinsVmVel{1}(:,1),saveBinsVmVel{1}(:,2))
    xline(saveBinsVmVel{1}(saveBinsVmVel{1}(:,2) == max(saveBinsVmVel{1}(:,2)),1),'-',num2str(saveBinsVmVel{1}(saveBinsVmVel{1}(:,2) == max(saveBinsVmVel{1}(:,2)),1)))
    legend('hide')
    title('T1 smoothVm')
    subplot(3,3,5);
    plot(saveBinsVmVel{2}(:,1),saveBinsVmVel{2}(:,2))
    xline(saveBinsVmVel{2}(saveBinsVmVel{2}(:,2) == max(saveBinsVmVel{2}(:,2)),1),'-',num2str(saveBinsVmVel{2}(saveBinsVmVel{2}(:,2) == max(saveBinsVmVel{2}(:,2)),1)))
    legend('hide')
    title('T2 smoothVm')
    subplot(3,3,8);
    plot(saveBinsVmVel{3}(:,1),saveBinsVmVel{3}(:,2))
    xline(saveBinsVmVel{3}(saveBinsVmVel{3}(:,2) == max(saveBinsVmVel{3}(:,2)),1),'-',num2str(saveBinsVmVel{3}(saveBinsVmVel{3}(:,2) == max(saveBinsVmVel{3}(:,2)),1)))
    legend('hide')
    title('T3 smoothVm')
    
    subplot(3,3,3);
    plot(saveBinsFRVel{1}(:,1),saveBinsFRVel{1}(:,2))
    xline(saveBinsFRVel{1}(saveBinsFRVel{1}(:,2) == max(saveBinsFRVel{1}(:,2)),1),'-',num2str(saveBinsFRVel{1}(saveBinsFRVel{1}(:,2) == max(saveBinsFRVel{1}(:,2)),1)))
    legend('hide')
    title('T1 FR')
    subplot(3,3,6);
    plot(saveBinsFRVel{2}(:,1),saveBinsFRVel{2}(:,2))
    xline(saveBinsFRVel{2}(saveBinsFRVel{2}(:,2) == max(saveBinsFRVel{2}(:,2)),1),'-',num2str(saveBinsFRVel{2}(saveBinsFRVel{2}(:,2) == max(saveBinsFRVel{2}(:,2)),1)))
    legend('hide')
    title('T2 FR')
    subplot(3,3,9);
    plot(saveBinsFRVel{3}(:,1),saveBinsFRVel{3}(:,2))
    xline(saveBinsFRVel{3}(saveBinsFRVel{3}(:,2) == max(saveBinsFRVel{3}(:,2)),1),'-',num2str(saveBinsFRVel{3}(saveBinsFRVel{3}(:,2) == max(saveBinsFRVel{3}(:,2)),1)))
    legend('hide')
    title('T3 FR')
    
    sgtitle(name)

    
% plots to check data quality of each trial
figure(18); clf;
f(1) = subplot(3,1,1);
    yyaxis left
    plot(processed_behaviourData{1}.time,processed_trialData{1}.smooth_Vm, 'k') 
    ylabel('Vm')
    hold on 
    yyaxis right
    plot(processed_behaviourData{1}.time,processed_behaviourData{1}.vel_dot, 'r')
f(2) = subplot(3,1,2);
    yyaxis left
    plot(processed_behaviourData{2}.time,processed_trialData{2}.smooth_Vm, 'k') 
    ylabel('Vm')
    hold on 
    yyaxis right
    plot(processed_behaviourData{2}.time,processed_behaviourData{2}.vel_dot, 'r')
f(3) = subplot(3,1,3);
    yyaxis left
    plot(processed_behaviourData{3}.time,processed_trialData{3}.smooth_Vm, 'k') 
    ylabel('Vm')
    hold on 
    yyaxis right
    plot(processed_behaviourData{3}.time,processed_behaviourData{3}.vel_dot, 'r')

figure(); clf;
f(4) = subplot(3,1,1);
    yyaxis left
    plot(processed_behaviourData{1}.time,processed_trialData{1}.scaledOutput_down, 'k') 
    ylabel('Vm')
    hold on 
    yyaxis right
    plot(processed_behaviourData{1}.time,processed_behaviourData{1}.vel_dot, 'r')
f(5) = subplot(3,1,2);
    yyaxis left
    plot(processed_behaviourData{2}.time,processed_trialData{2}.scaledOutput_down, 'k') 
    ylabel('Vm')
    hold on 
    yyaxis right
    plot(processed_behaviourData{2}.time,processed_behaviourData{2}.vel_dot, 'r')
f(6) = subplot(3,1,3);
    yyaxis left
    plot(processed_behaviourData{3}.time,processed_trialData{3}.scaledOutput_down, 'k') 
    ylabel('Vm')
    hold on 
    yyaxis right
    plot(processed_behaviourData{3}.time,processed_behaviourData{3}.vel_dot, 'r')
    
    linkaxes(f,'x');


%%
keep = input('Save? ','s');

corrData.corrSmoothVm = saveBinsVmVel;
corrData.corrFR = saveBinsFRVel;
corrData.corrScaledOutput = saveBinsSOVel;

    if strcmp(keep, 'y')
        cd(fileName) 
        savefig(h,'CrossCorr.fig')
        save('corrData.mat', 'corrData')
        cd('/Users/elenawesteinde/Documents/EphysCode') %'/Users/elenawesteinde/Documents/EphysCode''C:\Code\EphysCode'
    end
