function [rawData, trialData, trialMeta, behaviorData] = acquireOLTrial(trialRepeats, trialLength, inR, pattern_number, function_number, start_position)
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

%% Initialize settings
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
    if pattern_number == 1
        openLoop(settings.panels.barPattern, function_number, start_position);
    elseif pattern_number == 5
        openLoop(settings.panels.barPattern_prefHead, function_number, start_position);
    elseif pattern_number == 3
        openLoop(settings.panels.jumpPattern, function_number, start_position);
    end
    Panel_com('start');
    
    % Start Data Acquisition
    [rawDataTrial, trialTime] = niIO.startForeground;
    rawData{t} = rawDataTrial;
 
    % Data Processing
    [trialMeta.gain,trialMeta.mode,trialMeta.freq]= decodeTelegraphedOutput(rawDataTrial);
    
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
    h(3) = subplot(8,1,5);
    plot(behaviorData{t}.angle, 'k')
    ylabel('Pattern Angle')

    h(4) = subplot(8,1,6);
    plot( behaviorData{t}.vel_for, 'k')
    ylabel('Vel For')
    
    h(5) = subplot(8,1,7);
    plot( behaviorData{t}.vel_yaw, 'k')
    ylabel('Vel Yaw')
    
    h(6) = subplot(8,1,8);
    plot( behaviorData{t}.vel_side, 'k')
    ylabel('Vel Side')
    xlabel('Time (s)')

    linkaxes(h,'x');
    
    %plots ave Vm recorded when bar @ center of each panel
    if ~(function_number == 5)
        figure;
        medianfilteredOutput = medfilt1(trialData{t}.scaledOutput, 1000);
        framesX = behaviorData{t}.frX;
        edges = [0:8:96];
        [centers, mean_bin] = create_binned_mean(framesX, medianfilteredOutput, edges);
        plot(centers, mean_bin,'-o');
        xlabel('Frame number');
        ylabel('Vm');
    end
    
    
    Panel_com('stop');
    Panel_com('all_off');

end

trialMeta.trials      =  trialRepeats;
trialMeta.daqRate     =  niIO.Rate;
trialMeta.daqChIDs    = {niIO.Channels(:).ID};
trialMeta.daqChNames  = {niIO.Channels(:).Name};
trialMeta.pattern     =  pattern_number;

fprintf('\n******** acquireSimpleTrial Complete *********\n' )
end
 
function openLoop(pattern, func, start_position)
%% begins closedLoop setting in panels
    freq = 10;
    Panel_com('stop');
    %set pattern number
    Panel_com('set_pattern_id', pattern);
    %set open loop for x
    if pattern == 1
        Panel_com('set_mode', [4, 0]); 
        Panel_com('set_funcx_freq' , freq);
        Panel_com('set_posFunc_id', [1, func]);
        Panel_com('set_posFunc_id', [2, 0]);
        Panel_com('set_position', [start_position, 1]);
    elseif pattern == 5
        Panel_com('set_mode', [0, 0]); %[0,4]
        Panel_com('set_funcy_freq' , freq);
        Panel_com('set_posFunc_id', [2, func]);
        Panel_com('set_posFunc_id', [1, 0]);
        Panel_com('set_position', [start_position, 1]);
    elseif pattern == 3
        freq = 5; 
        Panel_com('set_mode', [0, 4]); %[0,4]
        Panel_com('set_funcy_freq' , freq);
        Panel_com('set_posFunc_id', [2, func]);
        Panel_com('set_posFunc_id', [1, 0]);
        Panel_com('set_position', [start_position, 1]);
    end

    %quiet mode on
    Panel_com('quiet_mode_on');
end

% functions should start from 0 and go from 0-95?