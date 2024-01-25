%% ephysSettings
% Specifies hard-coded parameters for ephys acquisition including data
% directory path and axopatch parameters. 
% EW 6/29/2020 JL 9/7/20

%% DEVICE CONFIGURATION
settings=struct();

% Set computer
comptype=computer;

% Set file directory, change path to your preference
settings.mainDataDir='E:\Dropbox (HMS)\Wilson_Lab_Data\ephys';

% NI Card device ID
settings.devID = 'Dev2';

% Set sampling rate
settings.sampRate=2e4; %Hz
settings.fictracRate = 60; %Hz

sensor_settings.sampRate = 4000;
sensor_settings.sensorPollFreq = 50;


%% NI BREAKOUT BOARD CHANNEL ASSIGNMENTS

% Break out box channels, change based on which channels on the DAQ
% breakout board you've plugged your BNC cables into
settings.bob.sinkCh = 0; 
settings.bob.scalCh = 8;
settings.bob.voltCh = 2;
settings.bob.currCh = 3;
settings.bob.gainCh = 4;
settings.bob.modeCh = 5;
settings.bob.freqCh = 6;

% Don't need to change this, normally
settings.bob.channelList = [settings.bob.sinkCh, settings.bob.scalCh, settings.bob.voltCh, settings.bob.currCh, settings.bob.gainCh, settings.bob.modeCh, settings.bob.freqCh]; 

% Break out box channels for behavior, changed based on which channels on
% the DAQ breakout board you've plugged your BNC cables into
settings.behavior.yawCh = 10;
settings.behavior.xCh = 11;
settings.behavior.yCh = 12;
settings.behavior.yawGainCh = 13;
settings.behavior.panelsXCh = 14;

settings.behavior.channelList = [settings.behavior.yawCh, settings.behavior.xCh, settings.behavior.yCh, settings.behavior.yawGainCh, settings.behavior.panelsXCh];


%% FORMAT OF THE RAW DATA. For use in DAQ settings.
% Don't need to change unless your change settings.bob.channelList or
% settings.behavior.channelList
settings.raw.sink = 1;
settings.raw.scaledOutput = 2;
settings.raw.voltage = 3;
settings.raw.current = 4;
settings.raw.gain = 5;
settings.raw.mode = 6;
settings.raw.frequency = 7;
settings.raw.yaw = 8;
settings.raw.x = 9;
settings.raw.y = 10;
settings.raw.yawgain = 11;
settings.raw.panels = 12;

settings.bob.aiType='SingleEnded'; % As opposed to 'differential' on the BOB - keep singleEnded.


%% PANELS SETTINGS
settings.panels.darkPattern = 1;
settings.panels.barPattern = 2;
settings.panels.numFrames = 96;

%% CURRENT AND VOLTAGE SETTINGS

% Current input settings - No signalCond at the moment
settings.current.betaRear  = 1; % Rear switch for current output set to beta = 1 mV/pA
settings.current.betaFront  = 1; % Front switch (CONFIG) for current output set to beta = 1      mV/pA
% settings.current.sigCond.Ch = 1;
% settings.current.sigCond.gain = 1;
% settings.current.sigCond.freq = 5; %kHz
MiliVOLTS_PER_VOLT = 1000; % 1000 mV/V  
settings.current.softGain   = MiliVOLTS_PER_VOLT/(settings.current.betaRear * settings.current.betaFront); % converted into pA and mV since 1pA/mV

% Voltage input settings - I am not using the signal conditioner currently
%settings.voltage.sigCond.Ch = 2;
%settings.voltage.sigCond.gain = 1;
%settings.voltage.sigCond.freq = 5; %kHz
settings.voltage.amplifierGain = 10; % 10 Vm, set coming out of the back of the amplifier
MiliVOLTS_PER_VOLT = 1000; % 1000 mV/V  
settings.voltage.softGain = MiliVOLTS_PER_VOLT / ( settings.voltage.amplifierGain); % To get voltage in mV

% Digital Voltage output settings:
%settings.daq.voltageDividerScaling = 0.0598; % voltage divider conversion factor, voltage divider cuts the volate by a factor of 0.0598
settings.daq.voltageDividerScaling = 1; % voltage divider removed from back of amplifier on 11/2
settings.daq.currentConversionFactor = 1 / (2000 * settings.current.betaFront * settings.daq.voltageDividerScaling); % V/pA   1 volt goes to 2 nA aka 2000 pA  
settings.daq.frontExtScale = 20 / 1000; %20mV/ 1000mV (1V) amplifier cuts the voltage down by this factor, every 1volt from the DAQ is 20mV into the Axopatch
settings.daq.voltageConversionFactor =  1 / (settings.daq.frontExtScale * settings.daq.voltageDividerScaling * MiliVOLTS_PER_VOLT); % use this for votlage clamp experiment commands 
%1 Volt = 2nA * Beta (1 normally)

% voltage clamp mode -
%   "20 mV/V", +1 V input produces +20 mV
%   input voltage range of -10 to +10 V produces -1000 to +1000 mV
settings.axopatch_mV_per_volt=20/1; % 20 mV / 1 V 

% current clamp mode -
%    "2 / (beta nA/V)" & beta = 1 nA/V
%    +1 V produces 2 nA (2000 pA)
%    input voltage range of -10 to +10 V produces -20 to +20 nA (-20,000 to +20,000 pA) 
settings.axopatch_picoAmps_per_volt=2000/1;

% If want to use -10 to +10 V (AO command) to inject -100 to +100 pA of current, 
% we need to scale the output of the daq board by a factor of 100/20,000=0.0050.
% voltage divider configuration is here:
% http://www.falstad.com/circuit/circuitjs.html?cct=$+1+0.000005+10.20027730826997+63+10+62%0Ar+288+336+288+256+0+4700%0Aw+288+256+368+160+0%0Aw+288+256+288+160+0%0Aw+288+96+288+64+0%0Aw+160+336+160+272+0%0Ag+432+96+464+64+0%0Aw+288+160+288+96+0%0AR+160+272+160+224+0+0+40+1+0+0+0.5%0Aw+160+336+160+400+0%0Aw+160+400+288+400+0%0Ar+288+400+288+336+0+4700%0Ar+368+160+400+128+0+22%0Ar+400+128+432+96+0+22%0A
settings.AO_output_scaling_factor=0.00476;
