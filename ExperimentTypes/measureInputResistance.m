function [inputData,inputResistance] = measureInputResistance()
% Aquires a trial of current clamp data obtained, with the extension command 
% on, and then use this response to calculate the input resistance. The 
% obtained trace and calculated input resistance are then saved

% INPUT
    % sampRate
 
% OUTPUT
    % inputData - obtained trace
    % inputResistance - calculated input resistance

%% Initialize Settings
ephysSettings
DAQsettings
% Prepare analog output channel
aO = niIO.addAnalogOutputChannel(devID,'ao0', 'Voltage');
aO.Name = 'External command';

%%
fprintf('\n^^^^^^^^^ Acquiring Input Resistance ^^^^^^^^^\n' )

% Set analog output command
currentPulse = [zeros(round(3*niIO.Rate),1); ones(round(3*niIO.Rate),1); zeros(round(1*niIO.Rate),1)]; % initialize
currentPulse = currentPulse * -5; % convert stim protocol to pA
currentPulse = currentPulse * (1/settings.axopatch_picoAmps_per_volt); % convert pA to Vout
currentPulse = currentPulse * (1/settings.AO_output_scaling_factor); % convert Vout to AO delivery from nidaq

output = zeros(size(currentPulse,1),1); % initialize
output(:,1) = currentPulse;

%% ACQUIRE TRIAL, CALCULATE INPUT RESISTANCE

niIO.queueOutputData(output);
rawData = niIO.startForeground;
rawData(:,1) = []; % remove sink channel

[trialMeta.gain,trialMeta.mode,trialMeta.freq]= decodeTelegraphedOutput(rawData);

% Process non-scaled data, change raw data channels based on your setup
inputData(:,1) = settings.current.softGain .* rawData(:,2);
inputData(:,2) = settings.voltage.softGain .* rawData(:,1);

switch trialMeta.mode
    % Voltage Clamp
    case {'Track','V-Clamp'}
        settings.scaledOutput.softGain = 1000 / (trialMeta.gain * settings.current.betaFront);
        inputData(:,3) = settings.scaledOutput.softGain .* rawData(:,4);  %mV
        % Current Clamp
    case {'I=0','I-Clamp Normal','I-Clamp Fast'}
        settings.scaledOutput.softGain = 1000 / (trialMeta.gain);
        inputData(:,3) = settings.scaledOutput.softGain .* rawData(:,4);  %pA
end

iDelta = mean(inputData((4*niIO.Rate):(6*niIO.Rate),1)) - mean(inputData(1:(2*niIO.Rate+1),1));
vDelta = mean(inputData((4*niIO.Rate):(6*niIO.Rate),3)) - mean(inputData(1:(2*niIO.Rate+1),3));
inputResistance = (vDelta*(10^-3))/(iDelta*(10^-12));

fprintf(['\n^^^^^^^^^^^^^^ Rinput = ' ,num2str(round(inputResistance/(10^6))),' MOhm ^^^^^^^^^^^^^\n'])
