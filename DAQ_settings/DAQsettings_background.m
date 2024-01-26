%% GENERATE DAQ SESSION AND ESTABLISH INPUT/OUTPUT CHANNELS

fprintf('\n************* Setting Parameters *************\n' )

ephysSettings;

niIO = daq.createSession('ni');
devID = settings.devID;

niIO.Rate = settings.sampRate;              % Sampling rate in Hz
niIO.DurationInSeconds = trialLength;       % Trial length in sec.

% Prepare analog input channels
AI = niIO.addAnalogInputChannel(devID, [settings.bob.channelList settings.behavior.channelList],'Voltage');
ai{settings.raw.sink} = 'Sink';
ai{settings.raw.gain} = 'Gain';   
ai{settings.raw.mode} = 'Mode';            % Amplifier mode
ai{settings.raw.frequency} = 'Frequency';       % Filter setting
ai{settings.raw.scaledOutput} = 'Scaled output';     
ai{settings.raw.voltage} = '10xVm';            % 10x fixed gain voltage output, no filter
ai{settings.raw.current} = 'Current';          % Current output filtered at 10kHz ?mV/pA
ai{settings.raw.yaw} = 'yaw';
ai{settings.raw.x} = 'x';
ai{settings.raw.y} = 'y';
ai{settings.raw.yawgain} = 'yaw gain';
ai{settings.raw.panelsX} = 'panels x';
ai{settings.raw.panelsY} = 'panels y';

for iAI = 1:length(ai)     % Assign amplifier analog input names
    AI(iAI).Name = ai{iAI};
    AI(iAI).TerminalConfig = settings.bob.aiType;
end

% Prepare analog output channel
aO = niIO.addAnalogOutputChannel(devID,'ao0', 'Voltage');
aO.Name = 'External command';

% To make handshaking robust use background acquisition with a logfile (read in later)
niIO.temp_log_filepath = fullfile(settings.mainDataDir,'logfile.dat');
niIO.temp_log_file_id = fopen(temp_log_filepath,'w+');
niIO.addlistener('DataAvailable',@(src,evt)utils.logDaqData(src, evt, temp_log_file_id, daq_exp_in_to_keep));


