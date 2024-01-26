function [rawData, trialData, trialMeta, behaviorData,dirname] = acquireCLTrial(trialRepeats, trialLength, inR)
% Acquisition code for baseline recording, no stimulus provided

% INPUT
    % trialLength - duration of each trial
% OUTPUT
    % trialData
    % trialMeta
    % rawData
    % behaviorData

%% Initialize settings
ephysSettings
[niIO, dataLog] = initalizeDAQ();
trial_complete = false;
exit_button_state = false; 

%% ACQUIRE TRIAL, SAVE EXPERIMENTAL PARAMETERS

% close fictrac if it's already running
closeFicTrac()

Panel_com('stop')

button = makeTerminateTrial();

% Start FicTrac in background from current experiment directory (config file must be in directory)
config = 'C:\Code\FicTrac\sample\config_current.txt'; 
dirname = 'C:\Users\ewest\OneDrive\Documents\Data\temp_fictrac';
foldername = fullfile(dirname,'fictracData');
mkdir(foldername);
configNew = fullfile(foldername,"config.txt");
copyfile(config, configNew);

expDir = foldername;
FT_PATH = 'C:\Code\FicTrac\bin\Release\fictrac.exe';
cmdStr = ['cd "', expDir, '" & start "" "',  FT_PATH, ...
        '" config.txt & exit'];
system(cmdStr);

%Panels Code
start = 1;
closedLoop(settings.panels.barPattern,start);
Panel_com('start');

trialMeta.pattern = settings.panels.barPattern;

% Call socket_client_CL_OL to open socket connection from fictrac to Phiget22 device
Socket_PATH = 'C:\Code\FicTrac\scripts';
SOCKET_SCRIPT_NAME = 'socket_client_360.py';
cmdstring = ['cd "' Socket_PATH '" & python ' SOCKET_SCRIPT_NAME ' &'];
[status] = system(cmdstring, '-echo');
    
       
% Start Data Acquisition

                fprintf(['\n************** Acquiring Trial *************\n'])

    startRecordingToLogFile(niIO)
    tic
    
    while ~exit_button_state && ~trial_complete
      % this loop will run until either the button GUI is pressed or the
      % max trial length has been reached
      % add a small delay so it doesn't try to run freakishly many loops
        pause(0.25)
        exit_button_state = button.EndTrialButton.Value; 
        if exist('trialLength', 'var')
            if trialLength <= toc
                trial_complete = true;
            end
        end
    end
    
    % close acquisition and save files
    stopRecordingToLogFile(niIO)
    dataLog = loadDaqDataFromLogFile(dataLog);
    rawData = dataLog.daq_data;
    clearLogFileData(dataLog);
    
    % close button 
    delete(button.UIFigure)
    clear button 
 
    % Data Processing
    [trialMeta.gain,trialMeta.mode,trialMeta.freq]= decodeTelegraphedOutput(rawData);
    
    % Process non-scaled data, adjust rawData channel based on which of
    % your channels are V vs I vs ScO
    current = settings.current.softGain .* rawData(:,settings.raw.current); % pA
    voltage = settings.voltage.softGain .* rawData(:,settings.raw.voltage); % mV
    
    time = (0:1/(settings.sampRate):(length(rawData(:,1))-1)/settings.sampRate)';
    
    trialData = table(time, current, voltage);
     
    switch trialMeta.mode
        % Voltage Clamp
        case {'Track','V-Clamp'}
            settings.scaledOutput.softGain = 1000 / (trialMeta.gain * settings.current.betaFront);
            trialData.scaledOutput = settings.scaledOutput.softGain .* adjustOffsetBasedGain(rawData(:,settings.raw.scaledOutput), trialMeta.gain);  %pA
            
            % Plot vclamp trial
            figure(99); clf;
            h(1) = subplot(4,1,1:3);
            plot(trialData.time, trialData.current, 'k')
            ylabel('Current (pA)')
            
            h(2) = subplot(4,1,4);
            plot(trialData.time, trialData.voltage, 'k')
            ylabel('Voltage (mV)')
            xlabel('Time (s)')
            
            linkaxes(h,'x')
            
        % Current Clamp
        case {'I=0','I-Clamp Normal','I-Clamp Fast'}
            settings.scaledOutput.softGain = 1000 / (trialMeta.gain);
            trialData.scaledOutput = settings.scaledOutput.softGain .* adjustOffsetBasedGain(rawData(:,settings.raw.scaledOutput), trialMeta.gain);  %mV
            
            % Plot iclamp trial
            figure(99); clf;
            h(1) = subplot(8,1,1);
            plot(trialData.time, trialData.current, 'k')
            ylabel('Current (pA)')
            
            h(2) = subplot(8,1,2:4);
            plot(trialData.time, trialData.scaledOutput, 'k')
            ylabel('Voltage (mV)')
            xlabel('Time (s)')
            
            sgtitle('CL_OL exp')
            hold on 
           
    end
    
    % Process behavior data
    [disp_for, disp_side, disp_yaw, frX, frY, angle, vel_for, vel_side, vel_yaw] = process_data(rawData, settings.panels.barPattern);
%     disp_for = disp_for';
%     disp_side = disp_side';
%     disp_yaw = disp_yaw';
%     vel_for = vel_for';
%     vel_side = vel_side';
%     vel_yaw = vel_yaw';

    behaviorData = table(disp_for, disp_side, disp_yaw, frX, frY, angle, vel_for, vel_side, vel_yaw); 
    
 time = (0:1/(settings.fictracRate/2):max(trialData.time))';
 %messy way to deal with unequal time vectors due to incomplete data
 %processing for these graphs, don't use this figure code for analysis
 
    % Plot behavior
    figure(99);
    h(3) = subplot(8,1,5);
    plot(trialData.time,behaviorData.angle, 'k')
    ylabel('Pattern Angle')

    h(4) = subplot(8,1,6);
    plot(trialData.time, behaviorData.vel_for, 'k')
    ylabel('Vel For (mm/s)')
    
    h(5) = subplot(8,1,7);
    plot(trialData.time, behaviorData.vel_yaw, 'k')
    ylabel('Vel Yaw (deg/s)')
    
    h(6) = subplot(8,1,8);
    plot(trialData.time, behaviorData.vel_side, 'k')
    ylabel('Vel Side (mm/s)')
    xlabel('Time (s)')

    linkaxes(h,'x');
    
%     figure(3);
%     %plot pattern fr vs Vm
%     medianfilteredOutput = medfilt1(trialData{t}.scaledOutput, 1000);
%     framesX = behaviorData{t}.frX;
%     edges = [0:8:96];
%     [centers, mean_bin] = create_binned_mean(framesX, medianfilteredOutput, edges);
%     plot(centers, mean_bin);
%     xlabel('Frame number');
%     ylabel('Vm');
    
    Panel_com('stop');
    Panel_com('all_off');
   
    % Kill FicTrac execution, sometimes it fails so keep doing it until it
    % shuts down
    
closeFicTrac()

trialMeta.daqRate     =  niIO.Rate;
trialMeta.daqChIDs    = {niIO.Channels(:).ID};
trialMeta.daqChNames  = {niIO.Channels(:).Name};

% Start FicTrac in background from current experiment directory (config file must be in directory)
runClosedLoopPattern()

fprintf('\n******** acquireCLTrial Complete *********\n' )
end