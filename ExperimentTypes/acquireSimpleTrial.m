function [rawData, trialData, trialMeta] = acquireSimpleTrial(trialRepeats, trialLength, inR)
% Acquisition code for baseline recording, no stimulus provided

% INPUT
    % trialRepeats - number of trial repeats
    % trialLength - duration of each trial
    % inR - [0/1] calculate input resistance?
% OUTPUT
    % trialData
    % trialMeta
    % rawData

%% Initialize Settings

daqreset;
ephysSettings;
DAQsettings;
% Prepare analog output channel
aO = niIO.addAnalogOutputChannel(devID,'ao0', 'Voltage');
aO.Name = 'External command';

%% CHECK INPUT RESISTANCE

if inR
    [~, trialMeta.inputR] = measureInputResistance();
else
    trialMeta.inputR = 'n/a';
end


%% ACQUIRE TRIAL, SAVE EXPERIMENTAL PARAMETERS

output = zeros(round(trialLength*niIO.Rate),1); % initialize

trialData = cell(1,trialRepeats);
rawData = cell(1, trialRepeats);

for t = 1:trialRepeats
    if trialRepeats>1
        fprintf(['\n************** Acquiring Trial ', num2str(t),' *************\n'])
    end

    niIO.queueOutputData(output);
    
    [rawDataTrial, trialTime] = niIO.startForeground;
    rawData{t} = rawDataTrial; 
    trialData{t}.time = trialTime;
    
    [trialMeta.gain, trialMeta.mode, trialMeta.freq]= decodeTelegraphedOutput(rawDataTrial);
    
    % Process non-scaled data, adjust rawData channel based on which of
    % your channels are V vs I vs ScO
    trialData{t}.current = settings.current.softGain .* rawDataTrial(:,settings.raw.current); % pA
    trialData{t}.voltage = settings.voltage.softGain .* rawDataTrial(:,settings.raw.voltage); % mV
    
    switch trialMeta.mode
        % Voltage Clamp
        case {'Track','V-Clamp'}
            settings.scaledOutput.softGain = 1000 / (trialMeta.gain * settings.current.betaFront);
            trialData{t}.scaledOutput = settings.scaledOutput.softGain .* rawDataTrial(:,settings.raw.scaledOutput);  %mV
            
            % Plot vclamp trial
            figure(86); clf;
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
            figure(86); clf;
            h(1) = subplot(4,1,1);
            plot(trialData{t}.time, trialData{t}.current, 'k')
            ylabel('Current (pA)')
            
            h(2) = subplot(4,1,2:4);
            plot(trialData{t}.time, trialData{t}.scaledOutput, 'k')
            ylabel('Voltage (mV)')
            xlabel('Time (s)')
            %ylim([-50 -30])
            
            sgtitle(['Trial ' num2str(t)])
            
            linkaxes(h,'x')
    end
end
trialMeta.trialDuration_s = length(trialData{t}.time)/niIO.Rate;
trialMeta.trials      =  trialRepeats;
trialMeta.daqRate     =  niIO.Rate;
trialMeta.daqChIDs    = {niIO.Channels(:).ID};
trialMeta.daqChNames  = {niIO.Channels(:).Name};

fprintf('\n******** acquireSimpleTrial Complete *********\n' )
 