function [rawData_all, trialData, trialMeta, behaviorData, dirname] = acquire_CLbehaviour_jump_Trial(trialRepeats, trialLength)
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
closedLoop_yFunc(settings.panels.barPattern, settings.panels.CL_jump_Function, start);

Panel_com('start');

trialMeta.pattern = settings.panels.barPattern; 
trialMeta.function = settings.panels.CL_jump_Function; 

% Call socket_client_360 to open socket connection from fictrac to Phiget22 device
Socket_PATH = 'C:\Code\FicTrac\scripts';
SOCKET_SCRIPT_NAME = 'socket_client_360.py';
cmdstring = ['cd "' Socket_PATH '" & python ' SOCKET_SCRIPT_NAME ' &'];
[status] = system(cmdstring, '-echo');

% Start Data Acquisition

trialData = cell(1,trialRepeats);
behaviorData = cell(1,trialRepeats);
rawData_all = cell(1,trialRepeats); 

for t = 1:trialRepeats
                fprintf(['\n************** Acquiring Trial *************\n'])
    button = makeTerminateTrial();      
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
    rawData_all{t} = rawData;
    clearLogFileData(dataLog);
    
    % close button 
    delete(button.UIFigure)
    clear button 
    
    time = (0:1/(settings.sampRate):(length(rawData(:,1))-1)/settings.sampRate)';
    
    
    
    % Process behavior data
    [disp_for, disp_side, disp_yaw, frX, frY, angle, vel_for, vel_side, vel_yaw] = process_data(rawData, settings.panels.barPattern);
    behaviorData{t} = table(disp_for, disp_side, disp_yaw, frX, frY, angle, vel_for, vel_side, vel_yaw); 
 
    % Plot behavior
    figure(99);
    h(1) = subplot(4,1,1);
    plot(time,behaviorData{t}.angle, 'k')
    ylabel('Pattern Angle')

    h(2) = subplot(4,1,2);
    plot(time, behaviorData{t}.vel_for, 'k')
    ylabel('Vel For (mm/s)')
    
    h(3) = subplot(4,1,3);
    plot(time, behaviorData{t}.vel_yaw, 'k')
    ylabel('Vel Yaw (deg/s)')
    
    h(4) = subplot(4,1,4);
    plot(time, behaviorData{t}.vel_side, 'k')
    ylabel('Vel Side (mm/s)')
    xlabel('Time (s)')

    linkaxes(h,'x');
    
    trial_complete = false;
    exit_button_state = false;
end  
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

% run CL pattern in between recorded trials
runClosedLoopPattern()

fprintf('\n******** acquireBehaviour_jump_Trial Complete *********\n' )
end
 