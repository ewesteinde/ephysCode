%% Load in exp data - use to reprocess all PFNd data acquired before 03/2021
close all 
clear

rootPath = '/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys'; %'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'; %change dep on comp
date = input('Date? ','s');
cell_num = input('Cell? ','s');
cell_num = strcat('cell_',cell_num);
trial = input('Trial? ','s');
trial = strcat('trial_',trial);
fileName = fullfile(rootPath,date,cell_num,trial);

cd(fileName)

load('trialData.mat');
load('trialMeta.mat');
load('behaviorData.mat');

if isfield(trialMeta, 'notes')
    disp(trialMeta.notes)
end

disp(trialMeta.trialType)

R = trialMeta.inputR/trialMeta.accessR;
disp(R)

if isfield(trialMeta, 'accessREnd')
    Rend = trialMeta.inputREnd/trialMeta.accessREnd;
    disp(Rend)
end

cd('/Users/elenawesteinde/Documents/EphysCode'); % 'C:\Code\JennyCode' %change dep on comp
%% Process Behaviour Data if good cell  

processed_behaviourData = cell(1,length(behaviorData));

for t = 1:length(behaviorData)
    ephysSettings
   %at this stage data has been lowpass filtered but is not downsampled or
   %in proper units
   %fictrac fps = 50 
 disp_for = behaviorData{t}.disp_for;
 disp_yaw = behaviorData{t}.disp_yaw;
 disp_side = behaviorData{t}.disp_side;
    
downsampled_disp_for = resample(disp_for,(25),settings.sampRate);
downsampled_disp_yaw = resample(disp_yaw,(25),settings.sampRate); 
downsampled_disp_side = resample(disp_side,(25),settings.sampRate); 

% smooth the position array
%filteredPosition = smoothdata(downsampled_cleanedPos,'rlowess',25);
% convert to proper units

filteredPosition_disp_for = ( downsampled_disp_for / (2*pi) ) * (pi*9);
filteredPosition_disp_yaw = ( downsampled_disp_yaw / (2*pi) ) * 360;
filteredPosition_disp_side = ( downsampled_disp_side / (2*pi) ) * (pi*9);

% differentiate into velocity
accumulatedPosition_for = filteredPosition_disp_for; 
accumulatedPosition_yaw = filteredPosition_disp_yaw; 
accumulatedPosition_side = filteredPosition_disp_side; 

velocity_for = gradient( filteredPosition_disp_for ) .* (25) ; 
velocity_yaw = gradient( filteredPosition_disp_yaw ) .* (25) ;
velocity_side = gradient( filteredPosition_disp_side ) .* (25) ;

% Calculate the distribution and take away values that are below 1% and above 99%
percentile25AV = prctile(velocity_for,1);
percentile975AV = prctile(velocity_for,99);
boundedVelocity_for = velocity_for;
boundedVelocity_for(velocity_for<percentile25AV | velocity_for>percentile975AV) = NaN;

percentile25AV = prctile(velocity_yaw,1);
percentile975AV = prctile(velocity_yaw,99);
boundedVelocity_yaw = velocity_yaw;
boundedVelocity_yaw(velocity_yaw<percentile25AV | velocity_yaw>percentile975AV) = NaN;

percentile25AV = prctile(velocity_side,1);
percentile975AV = prctile(velocity_side,99);
boundedVelocity_side = velocity_side;
boundedVelocity_side(velocity_side<percentile25AV | velocity_side>percentile975AV) = NaN;

% Linearly interpolate to replace the NaNs with values.
[pointsVectorAV] = find(~isnan(boundedVelocity_for));
valuesVectorAV = boundedVelocity_for(pointsVectorAV);
xiAV = 1:length(boundedVelocity_for);
interpVel_for = interp1(pointsVectorAV,valuesVectorAV,xiAV);

[pointsVectorAV] = find(~isnan(boundedVelocity_yaw));
valuesVectorAV = boundedVelocity_yaw(pointsVectorAV);
xiAV = 1:length(boundedVelocity_yaw);
interpVel_yaw = interp1(pointsVectorAV,valuesVectorAV,xiAV);

[pointsVectorAV] = find(~isnan(boundedVelocity_side));
valuesVectorAV = boundedVelocity_side(pointsVectorAV);
xiAV = 1:length(boundedVelocity_side);
interpVel_side = interp1(pointsVectorAV,valuesVectorAV,xiAV);

% smooth velocity
velocityOut_for = smoothdata(interpVel_for,'rlowess',15);
velocityOut_yaw = smoothdata(interpVel_yaw,'rlowess',15);
velocityOut_side = smoothdata(interpVel_side,'rlowess',15);

    %% Plot to check velocity & displacement outputs
    angle = resample(behaviorData{t}.angle,25,settings.sampRate);

    figure(1); clf;
        g(1) = subplot(4,1,1);
        plot( angle, 'k')
        ylabel('Pattern Angle')

        g(2) = subplot(4,1,2);
        plot( velocityOut_for, 'k')
        ylabel('Vel For')

        g(3) = subplot(4,1,3);
        plot( velocityOut_yaw, 'k')
        ylabel('Vel Yaw')

        g(4) = subplot(4,1,4);
        plot(velocityOut_side, 'k')
        ylabel('Vel Side')

        linkaxes(g,'x');

    figure(2); clf;
        g(1) = subplot(4,1,1);
        plot( angle, 'k')
        ylabel('Pattern Angle')

        g(2) = subplot(4,1,2);
        plot( accumulatedPosition_for, 'k')
        ylabel('disp For mm')

        g(3) = subplot(4,1,3);
        plot( accumulatedPosition_yaw, 'k')
        ylabel('disp Yaw deg')

        g(4) = subplot(4,1,4);
        plot( accumulatedPosition_side, 'k')
        ylabel('disp Side mm')

        linkaxes(g,'x');


    %% Flip v & disp side?
    fSide = input('flip side vel & disp? ','s');

    if strcmp(fSide, 'y')
        velocityOut_side = velocityOut_side.*-1;
        accumulatedPosition_side = accumulatedPosition_side.*-1;
% 
%         figure(3); clf;
%             g(1) = subplot(4,1,1);
%             plot(angle, 'k')
%             ylabel('Pattern Angle')
% 
%             g(2) = subplot(4,1,2);
%             plot(vel_for_down, 'k')
%             ylabel('Vel For')
% 
%             g(3) = subplot(4,1,3);
%             plot( vel_yaw_down, 'k')
%             ylabel('Vel Yaw')
% 
%             g(4) = subplot(4,1,4);
%             plot(vel_side_down, 'k')
%             ylabel('Vel Side')
% 
%             linkaxes(g,'x');
% 
%         figure(4); clf;
%             g(1) = subplot(4,1,1);
%             plot(behaviorData{t}.time, behaviorData{t}.angle, 'k')
%             ylabel('Pattern Angle')
% 
%             g(2) = subplot(4,1,2);
%             plot(behaviorData{t}.time, disp_for_mm, 'k')
%             ylabel('disp For')
% 
%             g(3) = subplot(4,1,3);
%             plot(behaviorData{t}.time, disp_yaw_deg, 'k')
%             ylabel('disp Yaw')
% 
%             g(4) = subplot(4,1,4);
%             plot(behaviorData{t}.time, disp_side_mm, 'k')
%             ylabel('disp Side')
%             xlabel('Time (s)')
% 
%             linkaxes(g,'x');
     end
%     
    
    %% Interpolate linearly so that original measurements are no longer stepwise, but has smooth transitions
    fHz = 1000;
    int_vel_for = resample(velocityOut_for,fHz,25);
    int_vel_yaw = resample(velocityOut_yaw,fHz,25);
    int_vel_side = resample(velocityOut_side,fHz,25);
    
    int_disp_for = resample(accumulatedPosition_for,fHz,25);
    int_disp_yaw = resample(accumulatedPosition_yaw,fHz,25);
    int_disp_side = resample(accumulatedPosition_side,fHz,25);
    
    angle = resample(behaviorData{t}.angle, 1000, settings.sampRate);
    angle = smoothdata(angle,'movmedian',20);
    angle(angle > 180) = 180;
    angle(angle < -180) = -180; 
    time = (0:1/fHz:max(behaviorData{t}.time))';
    
    

    figure(3); clf;
            g(1) = subplot(4,1,1);
            plot(time(1000:end-1000), angle(1000:end-1000), 'k')
            ylabel('Pattern Angle')

            g(2) = subplot(4,1,2);
            plot(time(1000:end-1000), int_vel_for(1000:end-1000), 'k')
            ylabel('Vel For mm/s')

            g(3) = subplot(4,1,3);
            plot(time(1000:end-1000), int_vel_yaw(1000:end-1000), 'k')
            ylabel('Vel Yaw deg/s')

            g(4) = subplot(4,1,4);
            plot(time(1000:end-1000), int_vel_side(1000:end-1000), 'k')
            ylabel('Vel Side mm/s')
            xlabel('Time (s)')
            
            linkaxes(g,'x');

    figure(4); clf;
            g(1) = subplot(4,1,1);
            plot(time(1000:end-1000), angle(1000:end-1000), 'k')
            ylabel('Pattern Angle')

            g(2) = subplot(4,1,2);
            plot(time(1000:end-1000), int_disp_for(1000:end-1000), 'k')
            ylabel('disp For rad/s')

            g(3) = subplot(4,1,3);
            plot(time(1000:end-1000), int_disp_yaw(1000:end-1000), 'k')
            ylabel('disp Yaw rad/s')

            g(4) = subplot(4,1,4);
            plot(time(1000:end-1000), int_disp_side(1000:end-1000), 'k')
            ylabel('disp Side rad/s')
            xlabel('Time (s)')
            
            linkaxes(g,'x');


    processed_behaviourData{t}.vel_for = int_vel_for(1000:end-1000);
    processed_behaviourData{t}.vel_yaw = int_vel_yaw(1000:end-1000);
    processed_behaviourData{t}.vel_side = int_vel_side(1000:end-1000);
    processed_behaviourData{t}.disp_for = int_disp_for(1000:end-1000);
    processed_behaviourData{t}.disp_yaw = int_disp_yaw(1000:end-1000);
    processed_behaviourData{t}.disp_side = int_disp_side(1000:end-1000);
    processed_behaviourData{t}.time = time(1000:end-1000);
    processed_behaviourData{t}.angle = angle(1000:end-1000);
    
    next = input('Next trial? ','s');
    if ~strcmp(next, 'y')
        error('fix issue & try again')
    end
end  

keep = input('Save? ','s');

if strcmp(keep, 'y')
    cd(fileName) 
    save('pro_behaviourData.mat','processed_behaviourData')
    cd('/Users/elenawesteinde/Documents/EphysCode')
end

%% Calculate firing rate
downsample_Hz = 1000; 
[fRate_sec, smoothVm, processed_trialData] =  calc_fRate(trialMeta,trialData, behaviorData, downsample_Hz, 6, 5, fileName);

%% Overlay raw Vm & Vf/Vs 

nActivity = processed_trialData{t}.scaledOutput_down;


figure(2);clf;      
        
        h(1) = subplot(3,1,1);
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.angle, 'k') 
        ylabel('angle')
       
        
        h(2) = subplot(3,1,2);
        yyaxis left
        plot(processed_behaviourData{t}.time,nActivity, 'k') 
        ylabel('Vm')
        %ylim([-55 -])
        hold on 
        yyaxis right
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_for, 'r')
        ylim([-5 5])
        ylabel('Vf mm/sec')
        
        h(3) = subplot(3,1,3);
        yyaxis left
        plot(processed_behaviourData{t}.time,nActivity, 'k')
        ylabel('Vm')
        %ylim([-55 -])
        hold on 
        yyaxis right
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_side, 'r')
        ylim([-5 5])
        ylabel('Vs mm/sec')
        
        
        linkaxes(h,'x');

       