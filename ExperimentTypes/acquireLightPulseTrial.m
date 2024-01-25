function [rawData, trialData, trialMeta] = acquireLightPulseTrial(trialRepeats, pulseLength, bufferLength, inR)
% Aquires a trial of current clamp data obtained, with the extension command 
% on, and then use this response to calculate the input resistance. The 
% obtained trace and calculated input resistance are then saved

% INPUT
    %amp, pulse amplitude
    %pulseLength stimulus duration
    %bufferLength time before and after stim before next trial
 
% OUTPUT
    % rawData
    % trialData
    % trialMeta
    
%% Initialize Settings
daqreset;
ephysSettings;
DAQsettings;
% Prepare analog output channel LIGHT PULSE
a1 = niIO.addAnalogOutputChannel(devID,'ao1', 'Voltage');
a1.Name = 'Shutter Pulse';

%% CHECK INPUT RESISTANCE

if inR
    [~, trialMeta.inputR] = measureInputResistance();

else
    trialMeta.inputR = 'n/a';
end

%% ACQUIRE TRIAL, SAVE EXPERIMENTAL PARAMETERS

trialData = cell(1,trialRepeats);

% Set analog output command
pulseLength = pulseLength/1000;   
lightPulse = [zeros(round(bufferLength*niIO.Rate),1); ones(round(pulseLength*niIO.Rate),1); zeros(round(bufferLength*niIO.Rate),1)]; % initialize
lightPulse = lightPulse * 5; % convert stim protocol to V

output = zeros(size(lightPulse,1),1); % initialize
output(:,1) = lightPulse;
trialMeta.trialAveBaseline = zeros(1,trialRepeats);
trialMeta.trialAveResp = zeros(1,trialRepeats);

for t = 1:trialRepeats
    if trialRepeats>1
        fprintf(['\n************** Acquiring Trial ', num2str(t),' *************\n'])
    end

    niIO.queueOutputData(output);

    % begin trial
    [rawData, trialTime] = niIO.startForeground;

    [trialMeta.gain, trialMeta.mode, trialMeta.freq]= decodeTelegraphedOutput(rawData);

    trialData{t}.time = trialTime;

    % Process non-scaled data, change rawData channels based on your setup
    trialData{t}.current = settings.current.softGain .* rawData(:,settings.raw.current); % pA
    trialData{t}.voltage = settings.voltage.softGain .* rawData(:,settings.raw.voltage); % mV

    trialData{t}.input = output(:,1)/5;

        switch trialMeta.mode
            % Voltage Clamp
            case {'Track','V-Clamp'}
                settings.scaledOutput.softGain = 1000 / (trialMeta.gain * settings.current.betaFront);
                trialData{t}.scaledOutput = settings.scaledOutput.softGain .* adjustOffsetBasedGain(rawData(:,settings.raw.scaledOutput), trialMeta.gain);  %pA
                trialMeta.trialAveBaseline(t) = sum(trialData{t}.scaledOutput(1:bufferLength*niIO.Rate))/length(trialData{t}.scaledOutput(1:bufferLength*niIO.Rate));
                trialMeta.trialResp(t) = -(max(abs(trialData{t}.scaledOutput(bufferLength*niIO.Rate:(bufferLength+.2)*niIO.Rate)))-abs(trialMeta.trialAveBaseline(t)));
            
                % Plot vclamp trial
                figure(11); clf;
                h(1) = subplot(5,1,1:3);
                plot(trialData{t}.time, trialData{t}.current, 'k')
                ylabel('Current (pA)')

                h(2) = subplot(5,1,4);
                plot(trialData{t}.time, trialData{t}.voltage, 'k')
                ylabel('Voltage (mV)')

                h(3) = subplot(5,1,5);
                plot(trialData{t}.time, trialData{t}.input,'k')
                ylabel('Light Stim')
                xlabel('Time (s)')
                
                linkaxes(h,'x')
                
                figure(22  ); clf;
                y = trialMeta.trialAveBaseline(1:t);
                yyaxis left
                plot(1:t,y,'-.')
                ylabel('Baseline (mV)')

                z = trialMeta.trialResp(1:t);
                yyaxis right
                plot(1:t,z,'-.')
                ylim([min(trialResp(1:t))-2 0])
                ylabel('Response (mV)')
                xlabel('Trial')

            % Current Clamp
            case {'I=0','I-Clamp Normal','I-Clamp Fast'}
                settings.scaledOutput.softGain = 1000 / (trialMeta.gain);
                trialData{t}.scaledOutput = settings.scaledOutput.softGain .* adjustOffsetBasedGain(rawData(:,settings.raw.scaledOutput), trialMeta.gain);  %pA  %pA
                trialMeta.trialAveBaseline(t) = sum(trialData{t}.scaledOutput(1:bufferLength*niIO.Rate))/length(trialData{t}.scaledOutput(1:bufferLength*niIO.Rate));
                trialMeta.trialResp(t) = -(max(abs(trialData{t}.scaledOutput(bufferLength*niIO.Rate:(bufferLength+.2)*niIO.Rate)))-abs(trialMeta.trialAveBaseline(t)));
                % Plot iclamp trial
                figure(11); clf;
                h(1) = subplot(5,1,1);
                plot(trialData{t}.time, trialData{t}.current, 'k');
                ylabel('Current (pA)')

                h(2) = subplot(5,1,2:4);
                plot(trialData{t}.time, trialData{t}.scaledOutput, 'k');
                ylabel('Voltage (mV)')

                h(3) = subplot(5,1,5);
                plot(trialData{t}.time, trialData{t}.input,'k');
                ylabel('Light Stim')
                xlabel('Time (s)')

                sgtitle(['Trial ' num2str(t)])
              
                linkaxes(h,'x')
                
                figure(22); clf; 
                y = trialMeta.trialAveBaseline(1:t);
                yyaxis left
                plot(1:t,y,'-.')
                hold on
                plot(1:t, movmean(trialMeta.trialAveBaseline(1:t),7),'k-');
                ylabel('Baseline (mV)')

                z = trialMeta.trialResp(1:t);
                
                yyaxis right
                plot(1:t,z,'-.')
                ylim([min(trialMeta.trialResp(1:t))-2 0])
                hold on
                plot(1:t, movmean(trialMeta.trialResp(1:t),7),'r-');
                ylabel('Response (mV)')
                xlabel('Trial')
        end

end

trialMeta.trialDuration_s = length(trialTime)/niIO.Rate; 
trialMeta.stimLength_ms = pulseLength;
trialMeta.trials      =  trialRepeats;
trialMeta.daqRate     =  niIO.Rate;
trialMeta.daqChIDs    = {niIO.Channels(:).ID};
trialMeta.daqChNames  = {niIO.Channels(:).Name};

fprintf('\n******** acquireLightPulseTrial Complete *********\n' )