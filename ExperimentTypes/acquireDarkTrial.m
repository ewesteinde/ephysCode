function [rawData, trialData, trialMeta, behaviorData] = acquireDarkTrial(trialRepeats, trialLength, inR)
% Acquisition code for baseline recording, no stimulus provided

% INPUT
    % trialRepeats - number of trial repeats
    % trialLength - duration of each trial
    % inR - [0/1] calculate input resistance?
% OUTPUT
    % trialData
    % trialMeta
    % rawData
    % behaviorData

% Initialize settings
daqreset;
ephysSettings;
DAQSettings_fictrac;
% Prepare analog output channel
aO = niIO.addAnalogOutputChannel(devID,'ao0', 'Voltage');
aO.Name = 'External command';

%% CHECK INPUT RESISTANCE

if inR
    [~, trialMeta.inputR] = measureInputResistance(niIO.Rate,gain);
else
    trialMeta.inputR = 'n/a';
end

%% ACQUIRE TRIAL, SAVE EXPERIMENTAL PARAMETERS

output = zeros(round(trialLength*niIO.Rate),1); % initialize
trialData = cell(1,trialRepeats);
behaviorData = cell(1, trialRepeats);
rawData = cell(1, trialRepeats);

for t = 1:trialRepeats
    if trialRepeats>1
        fprintf(['\n************** Acquiring Trial ', num2str(t),' *************\n'])
    end

    niIO.queueOutputData(output);
    
    % Panels code
    start = 1;
    closedLoop(settings.panels.darkPattern, start);
    Panel_com('start');
    
    % Start Data Acquisition
    [rawDataTrial, trialTime] = niIO.startForeground;
    rawData{t} = rawDataTrial;
 
    % Data Processing
    [trialMeta.gain, trialMeta.mode, trialMeta.freq]= decodeTelegraphedOutput(rawDataTrial);
    
    trialData{t}.time = trialTime;
    
    % Process non-scaled data, adjust rawData channel based on which of
    % your channels are V vs I vs ScO
    trialData{t}.current = settings.current.softGain .* rawDataTrial(:,settings.raw.current); % pA
    trialData{t}.voltage = settings.voltage.softGain .* rawDataTrial(:,settings.raw.voltage); % mV
    
    switch trialMeta.mode
        % Voltage Clamp
        case {'Track','V-Clamp'}
            settings.scaledOutput.softGain = 1000 / (trialMeta.gain * settings.current.betaFront);
            trialData{t}.scaledOutput = settings.scaledOutput.softGain .* adjustOffsetBasedGain(rawDataTrial(:,settings.raw.scaledOutput), trialMeta.gain);  %pA
            
            % Plot vclamp trial
            figure(1); clf;
            h(1) = subplot(4,1,1:3);
            plot(trialData{t}.time, trialData{t}.current, 'k')
            ylabel('Current (pA)')
            
            h(2) = subplot(4,1,4);
            plot(trialData{t}.time, trialData{t}.voltage, 'k')
            ylabel('Voltage (mV)')
            xlabel('Time (s)')
            
            linkaxes(h,'x')
            
        % Current Clamp
        case {'I=0','I-Clamp Normal','I-Clamp Fast'}
            settings.scaledOutput.softGain = 1000 / (trialMeta.gain);
            trialData{t}.scaledOutput = settings.scaledOutput.softGain .* adjustOffsetBasedGain(rawDataTrial(:,settings.raw.scaledOutput), trialMeta.gain);  %mV
            
            % Plot iclamp trial
            figure(1); clf;
            h(1) = subplot(8,1,1);
            plot(trialData{t}.time, trialData{t}.current, 'k')
            ylabel('Current (pA)')
            
            h(2) = subplot(8,1,2:4);
            plot(trialData{t}.time, trialData{t}.scaledOutput, 'k')
            ylabel('Voltage (mV)')
            xlabel('Time (s)')
            
            sgtitle(['Trial ' num2str(t)])
            hold on 
           
    end
    
    % Process behavior data
    pattern_num = 1;
    [disp_for, disp_side, disp_yaw, frX, frY, angle, vel_for, vel_side, vel_yaw] = process_data(rawDataTrial, pattern_num);
    behaviorData{t}.disp_for = disp_for;
    behaviorData{t}.disp_side = disp_side;
    behaviorData{t}.disp_yaw = disp_yaw;
    behaviorData{t}.frX = frX;
    behaviorData{t}.frY = frY;
    behaviorData{t}.angle = angle;
    behaviorData{t}.vel_for = vel_for;
    behaviorData{t}.vel_side = vel_side;
    behaviorData{t}.vel_yaw = vel_yaw;
    behaviorData{t}.time = trialTime;
    
    % Plot behavior
    time = (0:1/30:max(trialData{t}.time))';
    %messy way to deal with unequal time vectors due to incomplete data
    %processing for these graphs, don't use this figure code for analysis
 
    h(3) = subplot(8,1,5);
    plot(behaviorData{t}.time,behaviorData{t}.angle, 'k')
    ylabel('Pattern Angle')

    h(4) = subplot(8,1,6);
    plot( time, behaviorData{t}.vel_for, 'k')
    ylabel('Vel For (mm/s)')
    
    h(5) = subplot(8,1,7);
    plot( time, behaviorData{t}.vel_yaw, 'k')
    ylabel('Vel Yaw (deg/s)')
    
    h(6) = subplot(8,1,8);
    plot( time, behaviorData{t}.vel_side, 'k')
    ylabel('Vel Side (mm/s)')
    xlabel('Time (s)')

    linkaxes(h,'x');
    
    Panel_com('stop');
    Panel_com('all_off');

end

trialMeta.trials      =  trialRepeats;
trialMeta.daqRate     =  niIO.Rate;
trialMeta.daqChIDs    = {niIO.Channels(:).ID};
trialMeta.daqChNames  = {niIO.Channels(:).Name};

fprintf('\n******** acquireSimpleTrial Complete *********\n' )
end
 
function closedLoop(pattern, startPosition)
   %% begins closedLoop setting in panels
    Panel_com('stop');
    %set arena
    pause(.05)
    Panel_com('set_config_id', 1);
    %set brightness level
    pause(.05)
%     Panel_com('g_level_1');
    %set pattern number
    Panel_com('set_pattern_id', pattern);
    Panel_com('set_position', [96, startPosition]);
    %% 
    %set closed loop for x
    pause(.05)
    Panel_com('set_mode', [3, 0]);
    pause(.05)
    Panel_com('quiet_mode_on');
end