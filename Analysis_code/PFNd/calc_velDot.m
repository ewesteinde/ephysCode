%% calculate the preferred translation velocity of the PFNd neuron
% 
% clear
% rootPath = '/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys'; % ;'/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys' 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'%change dep on comp
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
% load('trialMeta.mat');
% load('trialData.mat');
% 
% if isfield(trialMeta, 'notes')
%     disp(trialMeta.notes)
% end
% 
% cd('/Users/elenawesteinde/Documents/EphysCode'); %  'C:\Code\EphysCode' '/Users/elenawesteinde/Documents/EphysCode'%change dep on comp

%%
function [processed_behaviourData] = calc_velDot(offset, processed_trialData, processed_behaviourData, fileName)

scaledOutput = [];
fRate_sec = [];
angle = [];
Vf = [];
Vs = [];
time = [];

for t = 1:length(processed_behaviourData)

%     if trialMeta.fly.timestamp < '26-Mar-2021' && trialMeta.fly.timestamp > '01-Jan-2021'
%         downsample_Hz = 1000; 
%         ephysSettings
%         activity = resample(trialData{t}.scaledOutput, downsample_Hz, settings.sampRate);
%         %correct for leak junction potential 
%         activity = activity-13;
%         scaledOutput = activity;
%     else
        scaledOutput = processed_trialData{t}.scaledOutput_down; 
%     end
    
    angle = processed_behaviourData{t}.angle;
    Vf =  processed_behaviourData{t}.vel_for;
    Vs = processed_behaviourData{t}.vel_side;

time = processed_behaviourData{t}.time;
Vf = Vf';
Vs = Vs';

Vs_nopref = Vs.*cosd(90 - offset);
Vf_nopref = Vf.*cosd(offset);
velDot = Vs_nopref + Vf_nopref;

processed_behaviourData{t}.vel_dot = velDot;
    
figure();clf;  
    h(1) = subplot(2,1,1);
    plot(time, angle, 'k') 
    ylabel('angle')


    h(2) = subplot(2,1,2);
    yyaxis left
    plot(time,scaledOutput, 'k') 
    ylabel('Vm')
    %ylim([-70 -30])
    hold on 
    yyaxis right
    plot(time,processed_behaviourData{t}.vel_dot, 'r')
    ylim([-15 15])

    ylabel('Pref translational vel mm/sec')
    
    linkaxes(h,'x');
 

    next = input('Next trial? ','s');
    if ~strcmp(next, 'y')
        error('fix issue & try again')
    end

end

keep = input('Save? ','s');

if strcmp(keep, 'y')
    cd(fileName) 
    save('pro_behaviourData.mat','processed_behaviourData')
    cd('/Users/elenawesteinde/Documents/EphysCode') %'/Users/elenawesteinde/Documents/EphysCode/Analysis_code' 'C:\Code\EphysCode'
end
end

% create heatmap to est heading pref of cell based pref vel vector
           
 %% calculate velocity of fly, dep on heading preference of cell

 
% Vspref = Vs.*cosd((angle + 90)-wrapTo180(prefHead+offset));
% Vfpref = Vf.*cosd(angle-wrapTo180(prefHead+offset)); 
%  
% 
% vel_noPref = Vs_nopref + Vf_nopref;
% 
% h(4) = subplot(4,1,4);
%     yyaxis left
%     plot(time,scaledOutput, 'k') 
%     ylabel('Vm')
%     ylim([-70 -30])
%     hold on 
%     yyaxis right
%     plot(time,vel_noPref, 'r')
%     ylim([-12 12])
%     ylabel('Speed mm/sec')
%     
%     linkaxes(h,'x');
%     

