%% load in experiment 

clear
close all

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
openfig('CrossCorr.fig');

if isfield(trialMeta, 'notes')
    disp(trialMeta.notes)
end


cd('/Users/elenawesteinde/Documents/EphysCode'); %  'C:\Code\EphysCode' '/Users/elenawesteinde/Documents/EphysCode'%change dep on comp

 %  'C:\Code\EphysCode' '/Users/elenawesteinde/Documents/EphysCode'%change dep on comp

% extract timepoints when the fly is moving
minVel = 0.5; %mm/s only taking into account side & forward velocities
trial_fragments = cell(1,3); 

for t = 1:length(processed_trialData)
    speed_singleTrial = abs(processed_behaviourData{t}.vel_dot); 
    [index] = find(speed_singleTrial > minVel); 
    count = 1;
    start = 0;
    frag = {};
    shiftPoint = [];
    for i = 2:length(index)
        if index(i) - index(i-1) > 20
            shiftPoint(1,2) = i-1;
            if count == 1
                frag{count} = index(1:shiftPoint(1,2)); 
                count = count+1;
                shiftPoint(1,1) = i;
            else
                frag{count} = index(shiftPoint(1,1):shiftPoint(1,2));
                count = count+1;
                shiftPoint(1,1) = i;
            end 
        end

        if i == length(index)
            if isempty(shiftPoint)
                frag{count} = index(1:i);
            else
                frag{count} = index(shiftPoint(1,1):i);
            end
        end
    end
    
%     for c = 1:length(frag)
%         if frag{c+1}(1) - frag{c}(end) < 500
%             cont_frag{c} = 
    
    trial_fragments{t} = frag; 
end

%

name = strcat(date,', ',cell_num,', ',trial);
      

ephysSettings;

count = 1; 

xcorrSum = cell(1,3); 
xcorrAve = cell(1,3); 

for t = 1:length(trial_fragments)
    velDot = processed_behaviourData{t}.vel_dot;  
%     smoothVm = processed_trialData{t}.smooth_Vm;
%     smVm = smoothdata(smoothVm,'loess',100);
    smVm = processed_trialData{t}.fRate_sec;
    nanIDX  = find( isnan(velDot) );
    xcorrSum{t} = zeros(1001,1);

    if ~isempty(nanIDX)
        if nanIDX(1) == 1
            velDot = velDot(nanIDX(end)+1:length(velDot)); 
            smVm = smVm(nanIDX(end)+1:length(smVm)); 
        elseif nanIDX(end) == length(velDot)
            velDot = velDot(1:nanIDX(1)-1);
            smVm = smVm(1:nanIDX(1)-1); 
        else
            error('NaN values inside velDot')
        end
    end
    count = 0; 
    timeUsed = 0; 
    
    for i = 1:length(trial_fragments{t})
        if length(trial_fragments{t}{i}) >= 1000
            timeUsed = timeUsed + length(trial_fragments{t}{i}); 
            velDot_frag = velDot(trial_fragments{t}{i}(1):trial_fragments{t}{i}(end));
            smVm_frag = smVm(trial_fragments{t}{i}(1):trial_fragments{t}{i}(end));

            [cVm_frag, lagsVm_frag] = xcorr(smVm_frag',velDot_frag,500);
% 
%             figure(2); clf; 
%             yyaxis left
%             plot(velDot_frag)
%             yyaxis right
%             plot(smVm_frag)

            xcorrSum{t} = xcorrSum{t} + cVm_frag; 

            count = count + 1;
        end
    end
    xcorrAve{t} = xcorrSum{t}/count;
%     figure(1);clf
%     plot(lagsVm_frag, xcorrAve{t})
     xcorrLag = lagsVm_frag; 
end
%
h = figure(20); clf; 


    subplot(3,1,1);
    plot(xcorrLag,xcorrAve{1})
    xline(xcorrLag(xcorrAve{1} == max(xcorrAve{1})),'-',num2str(xcorrLag(xcorrAve{1} == max(xcorrAve{1}))))
    legend('hide')
    title('T1 fRate')
    subplot(3,1,2);
    plot(xcorrLag,xcorrAve{2})
    xline(xcorrLag(xcorrAve{2} == max(xcorrAve{2})),'-',num2str(xcorrLag(xcorrAve{2} == max(xcorrAve{2}))))
    legend('hide')
    title('T2 fRate')
    subplot(3,1,3);
    plot(xcorrLag,xcorrAve{3})
    xline(xcorrLag(xcorrAve{3} == max(xcorrAve{3})),'-',num2str(xcorrLag(xcorrAve{3} == max(xcorrAve{3}))))
    legend('hide')
    title('T3 fRate')
    
    sgtitle(name)
    
    timeUsedsec = timeUsed/1000; 


keep = input('Save? ','s');

    if strcmp(keep, 'y')
        cd(fileName) 
        savefig(h,'CrossCorr_no0Vel.fig')
        cd('/Users/elenawesteinde/Documents/EphysCode') %'/Users/elenawesteinde/Documents/EphysCode''C:\Code\EphysCode'
    end

%% plots to check data quality of each trial
% figure(18); clf;
% f(1) = subplot(3,1,1);
%     yyaxis left
%     plot(processed_behaviourData{1}.time,processed_trialData{1}.smooth_Vm, 'k') 
%     ylabel('Vm')
%     hold on 
%     yyaxis right
%     plot(processed_behaviourData{1}.time,processed_behaviourData{1}.vel_dot, 'r')
% f(2) = subplot(3,1,2);
%     yyaxis left
%     plot(processed_behaviourData{2}.time,processed_trialData{2}.smooth_Vm, 'k') 
%     ylabel('Vm')
%     hold on 
%     yyaxis right
%     plot(processed_behaviourData{2}.time,processed_behaviourData{2}.vel_dot, 'r')
% f(3) = subplot(3,1,3);
%     yyaxis left
%     plot(processed_behaviourData{3}.time,processed_trialData{3}.smooth_Vm, 'k') 
%     ylabel('Vm')
%     hold on 
%     yyaxis right
%     plot(processed_behaviourData{3}.time,processed_behaviourData{3}.vel_dot, 'r')
% 
% figure(); clf;
% f(4) = subplot(3,1,1);
%     yyaxis left
%     plot(processed_behaviourData{1}.time,processed_trialData{1}.scaledOutput_down, 'k') 
%     ylabel('Vm')
%     hold on 
%     yyaxis right
%     plot(processed_behaviourData{1}.time,processed_behaviourData{1}.vel_dot, 'r')
% f(5) = subplot(3,1,2);
%     yyaxis left
%     plot(processed_behaviourData{2}.time,processed_trialData{2}.scaledOutput_down, 'k') 
%     ylabel('Vm')
%     hold on 
%     yyaxis right
%     plot(processed_behaviourData{2}.time,processed_behaviourData{2}.vel_dot, 'r')
% f(6) = subplot(3,1,3);
%     yyaxis left
%     plot(processed_behaviourData{3}.time,processed_trialData{3}.scaledOutput_down, 'k') 
%     ylabel('Vm')
%     hold on 
%     yyaxis right
%     plot(processed_behaviourData{3}.time,processed_behaviourData{3}.vel_dot, 'r')
%     
%     linkaxes(f,'x');


%%
