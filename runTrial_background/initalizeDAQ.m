% initilaize DAQ session function, compare to my code
%niIO should replace all calls of niIO
function [niIO, dataLog] = initalizeDAQ()    
    
    % Experimental Signals In
    daq.reset();
    ephysSettings
    niIO = daq.createSession('ni');
    devID = settings.devID;
    
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
    
    niIO.IsContinuous = true;
    
    % To make handshaking robust use background acquisition with a logfile (read in later)
    dataLog.temp_log_filepath = fullfile(settings.mainDataDir,'logfile.dat');
    dataLog.temp_log_file_id = fopen(dataLog.temp_log_filepath,'w+');
    % when the event 'DataAvailable' occurs I belive this will execute the logDaqData
    % function, need to figure out where daq_exp_in_to_keep came from in
    % Stephen's code
    niIO.addlistener('DataAvailable',@(src,evt)logDaqData(src,evt,dataLog.temp_log_file_id));
    % Set acquisition rate
    niIO.Rate = settings.sampRate;

end


