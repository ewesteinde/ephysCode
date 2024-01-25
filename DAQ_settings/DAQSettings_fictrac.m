  
%% GENERATE DAQ SESSION AND ESTABLISH INPUT/OUTPUT CHANNELS

fprintf('\n************* Setting Parameters *************\n' )


niIO = daq.createSession('ni');
devID = settings.devID;

niIO.Rate = settings.sampRate;              % Sampling rate in Hz

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


