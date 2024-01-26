
function [ velocityOut , accumulatedPositionOut ] = ficTracSignalDecoding_noSmooth( ficTracBallPosition , sampleRate , FicTracRate, yaw)
%FICTRACSIGNALDECODING takes a fictrac position value and extracts ball velocity 
%  This function will take a FicTrac output signal as aquired by the DAQ
%  as an analog signal and solve for the ball's angular velocity in the given 
%  dimention. To do this the signal is then UNWRAPPED to handle the abrupt
%  transitions caused by when the ball rotates completely and the signal
%  resets (0->10 volts or 10 volts -> 0) transistions.  Then the signal will be
%  further CLEANED to remove extra position values surrounding those signal
%  reset time points. Then the position signal is LOW PASS FILTERED to
%  remove noise. Then the velocity of the ball with be solved for in
%  degree/s by using the diff function, and taking into consideration the  
%  sample rate the data was collected at.  Velocity values above what a
%  resonable fly would turn the ball at (maxFlyVelocity) are discarded
%  
%   INPUTS
%   ficTracBallPosition  - array containing data from 0-10 volts relating
%   to the balls position
%
%   sampleRate - Rate the data was aquired at ( samples/ second )
%   
%   lowPassFilterCutOff - frequency that the position signal will be low
%   pass filtered at (Hz)
%
%   maxFlyVelocity - max value of realistic fly movement (deg/s) 
%
%   OUTPUT
%   velocityOut -array containing ball's instentanous velocity (degree/sec)
%   accumulatedPositionOut - array containing the filtered and unwraped
%   position signal
%
%   Yvette Fisher 1/2018
%   Modified by Jenny Lu 2/2019 for better derivative -- 3/2019 reverted to
%   butterworth filter for consistency
%   EW edits 3/23/21 to output Yaw rotational velocity
%   in degrees/sec. Forward & side in mm/s, changed smoothing steps
% ------------------------------
FICTRAC_MAX_VOLTAGE = 10;  % volts

% transfrom ficTrac signal into radians  
posRadians = ficTracBallPosition .* 2 .* pi ./ FICTRAC_MAX_VOLTAGE; 

% upwrap position signal
unwrappedPos = unwrap( posRadians );

% find indexes where the unwrapping happened (tolerace = pi)
upwrappedIndexes = find ( abs( diff( posRadians )) > pi); 

NUM_SAMPLES_FROM_WRAP_TO_REPLACE = 2;
% handle edge case so we don't fill off the edge of the trace
upwrappedIndexes = upwrappedIndexes( upwrappedIndexes > NUM_SAMPLES_FROM_WRAP_TO_REPLACE & upwrappedIndexes < (length ( unwrappedPos ) - NUM_SAMPLES_FROM_WRAP_TO_REPLACE) ); 

cleanedPos = unwrappedPos;
% replace potentially problematic indexes with Nan
for i = 1: length ( upwrappedIndexes )
    index_start = upwrappedIndexes(i) -  NUM_SAMPLES_FROM_WRAP_TO_REPLACE ; 
    index_end = upwrappedIndexes(i) +  NUM_SAMPLES_FROM_WRAP_TO_REPLACE ; 
    
    cleanedPos( index_start : index_end ) = NaN;
end

% replace NaN values with the last preceding value that was a real number
nanIDX = find( isnan( cleanedPos ) ); % find NaN indexes
% replace with preceeding value
while( ~isempty( nanIDX ) )
    cleanedPos(nanIDX) = cleanedPos(nanIDX - 1);
    
    % find any remaining NaN
    nanIDX  = find( isnan(cleanedPos) );
end

% downsample position data to match FicTrac's output
[downsampled_cleanedPos,~] = resample_new(cleanedPos,(FicTracRate/2),sampleRate); 

% smooth the position array
filteredPosition = downsampled_cleanedPos; 
%smoothdata(downsampled_cleanedPos,'rlowess',10);

% convert to proper units
 if yaw == 1
        filteredPosition_unitConv = ( filteredPosition / (2*pi) ) * 360;
    else
        filteredPosition_unitConv = ( filteredPosition / (2*pi) ) * (pi*9);
end


% differentiate into velocity
accumulatedPositionOut = filteredPosition_unitConv; 

velocity = gradient( filteredPosition ) .* (FicTracRate/2) ; 

 if yaw == 1
        velocity_unitConv = ( velocity / (2*pi) ) * 360;
    else
        velocity_unitConv = ( velocity / (2*pi) ) * (pi*9);
end
 

% Calculate the distribution and take away values that are below 2.5% and above 97.5%
percentile25AV = prctile(velocity_unitConv,1);
percentile975AV = prctile(velocity_unitConv,99);
boundedVelocity = velocity_unitConv;
boundedVelocity(velocity_unitConv<percentile25AV | velocity_unitConv>percentile975AV) = NaN;

% Linearly interpolate to replace the NaNs with values.
[pointsVectorAV] = find(~isnan(boundedVelocity));
valuesVectorAV = boundedVelocity(pointsVectorAV);
xiAV = 1:length(boundedVelocity);
interpVel = interp1(pointsVectorAV,valuesVectorAV,xiAV);

% smooth velocity
%velocityOut = smoothdata(velocity,'gaussian',10);
velocityOut = interpVel;

%% plotting to check how well unwrapping, cleaning and filtering worked
% can be commented out once you are happy with the parameters
 
figure('Position',[50 50 800 900]),
h(1) = subplot(6,1,1);
plot(ficTracBallPosition)
title('Raw voltage signal');
ylim([0 10]);
xlim([0 length(ficTracBallPosition)]);
set(gca,'xticklabel',{[]})

h(2) = subplot(6,1,2);
plot(posRadians)
title('Signal in radians');
ylim([0 2*pi]);
xlim([0 length(posRadians)]);
set(gca,'xticklabel',{[]})

h(3) = subplot(6,1,3);
plot(cleanedPos)
title('Unwrapped signal in radians');
xlim([0 length(cleanedPos)]);
set(gca,'xticklabel',{[]})

h(4) = subplot(6,1,4);
plot(downsampled_cleanedPos)
title('Downsampled signal in radians');
xlim([0 length(downsampled_cleanedPos)]);
set(gca,'xticklabel',{[]})

h(5) = subplot(6,1,5);
plot(filteredPosition_unitConv)
title('lowpass Filtered position signal in mm or deg');
xlim([0 length(filteredPosition_unitConv)]);
set(gca,'xticklabel',{[]})

h(6) = subplot(6,1,6);
plot(velocity)
title('velocity signal in mm or deg');
xlim([0 length(velocity)]);
set(gca,'xticklabel',{[]})


%linkaxes(h,'x');