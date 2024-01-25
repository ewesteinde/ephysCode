function [rawData,trialData, trialMeta] = acquireCurrentPulseTrial(trialRepeats, amp, pulseLength, bufferLength, inR)
% Aquires a trial of current clamp data obtained, with the extension command 
% on, and then use this response to calculate the input resistance. The 
% obtained trace and calculated input resistance are then saved

% INPUT
    % amp, pulse amplitude
    %pulseLength stimulus duration
    %bufferLength time before and after stim before next trial
 
% OUTPUT
    % inputData - obtained trace
    % inputResistance - calculated input resistance
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

trialData = cell(1,trialRepeats);
for t = 1:trialRepeats
    if trialRepeats>1
        fprintf(['\n************** Acquiring Trial ', num2str(t),' *************\n'])
    end
    
% Set analog output command
currentPulse = [zeros(round(bufferLength*niIO.Rate),1); ones(round(pulseLength*niIO.Rate),1); zeros(round(bufferLength*niIO.Rate),1)]; % initialize
currentPulse = currentPulse * amp; % convert stim protocol to pA
currentPulse = currentPulse * (1/settings.axopatch_picoAmps_per_volt); % convert pA to Vout
currentPulse = currentPulse * (1/settings.AO_output_scaling_factor); % convert Vout to AO delivery from nidaq

output = zeros(size(currentPulse,1),1); % initialize
output(:,1) = currentPulse;

niIO.queueOutputData(output);
rawData = niIO.startForeground;
rawData(:,1) = []; % remove sink channel

[trialMeta.gain,trialMeta.mode,trialMeta.freq]= decodeTelegraphedOutput(rawData);

trialTime = (1:1:length(output))/niIO.Rate;

% Process non-scaled data, change rawData channels based on your setup
trialData{t}(:,1) = settings.current.softGain .* rawData(:,2);
trialData{t}(:,2) = settings.voltage.softGain .* rawData(:,1);

    switch trialMeta.mode
        % Voltage Clamp
        case {'Track','V-Clamp'}
            settings.scaledOutput.softGain = 1000 / (trialMeta.gain * settings.current.betaFront);
            trialData{t}(:,3) = settings.scaledOutput.softGain .* rawData(:,4);  %mV
            
            % Plot vclamp trial
            figure(1); clf;
            h(1) = subplot(4,1,1:3);
            plot(trialTime', trialData{t}(:,3), 'k')
            ylabel('Current (pA)')
            
            h(2) = subplot(4,1,4);
            plot(trialTime', trialData{t}(:,2), 'k')
            ylabel('Voltage (mV)')
            xlabel('Time (s)')
            
            linkaxes(h,'x')
            
        % Current Clamp
        case {'I=0','I-Clamp Normal','I-Clamp Fast'}
            settings.scaledOutput.softGain = 1000 / (trialMeta.gain);
            trialData{t}(:,3) = settings.scaledOutput.softGain .* rawData(:,4);  %pA
            
            % Plot iclamp trial
            figure(1); clf;
            h(1) = subplot(4,1,1);
            plot(trialTime', trialData{t}(:,1), 'k')
            ylabel('Current (pA)')
            
            h(2) = subplot(4,1,2:4);
            plot(trialTime', trialData{t}(:,3), 'k')
            ylabel('Voltage (mV)')
            xlabel('Time (s)')
            
            sgtitle(['Trial ' num2str(t)])
            
            linkaxes(h,'x')
    end

%% troubleshooting, temporary code, gives change in I vs V traces during I pulse, verify stim amp is true

iDelta = mean(trialData{t}(((bufferLength+1)*niIO.Rate):((bufferLength+pulseLength)*niIO.Rate),1)) - mean(trialData{t}(1:((bufferLength-1)*niIO.Rate+1),1));
vDelta = mean(trialData{t}(((bufferLength+1)*niIO.Rate):((bufferLength+pulseLength)*niIO.Rate),3)) - mean(trialData{t}(1:((bufferLength-1)*niIO.Rate+1),3));
IV(t,1) = iDelta;
IV(t,2) = vDelta;

end

trialMeta.trialDuration_s = length(trialTime)/niIO.Rate; 
trialMeta.stimLength_s = pulseLength;
trialMeta.stimAmp_pA = amp;
trialMeta.trials      =  trialRepeats;
trialMeta.daqRate     =  niIO.Rate;
trialMeta.daqChIDs    = {niIO.Channels(:).ID};
trialMeta.daqChNames  = {niIO.Channels(:).Name};

fprintf('\n******** acquireCurrentPulseTrial Complete *********\n' )