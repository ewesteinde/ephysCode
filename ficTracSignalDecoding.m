% <<<<<<< HEAD
function [ velocityOut , accumulatedPositionOut ] = ficTracSignalDecoding( ficTracBallPosition , sampleRate , FicTracRate, yaw)
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
if yaw 
    for_pos = ficTracBallPosition(:,2);
    ficTracBallPosition = ficTracBallPosition(:,1);
end

% transfrom ficTrac signal into radians  
posRadians = ficTracBallPosition .* 2 .* pi ./ FICTRAC_MAX_VOLTAGE; 

% upwrap position signal
unwrappedPos = unwrap( posRadians );

% when fictrac breaks or before the fly moves it creates a partial wrap of the signal that unwrap
% doesn't fix, tends to happen in the yaw channel --> do it manually
if yaw
    setTo0 = 0;
    BadIndexes = find(abs(gradient(unwrappedPos)) > 0.2); 
    if ~isempty(BadIndexes)
        sigChange_idx = ischange(for_pos);
        startMovIdx = find(sigChange_idx,1,'first'); 
        if BadIndexes(1) < startMovIdx
            unwrappedPos(1:startMovIdx) = 0;
            setTo0 = 1;
        elseif isempty(startMovIdx) %fly didn't move for the entire trial 
             unwrappedPos(1:end) = 0; 
             setTo0 = 1; 
        end
        try
            prc99 = prctile(abs(gradient(unwrappedPos)),99.99);
            if prc99 > 0.1
                prc99 = 0.05; % rough catch for when there was fictrac issues & fly didn't move much
            end
            BadIndexes = find(abs(gradient(unwrappedPos)) > 0.2); 
            BadIndexes = find(abs(gradient(unwrappedPos)) > prc99); 
            fixedPos = unwrappedPos;
            NUM_SAMPLES_FROM_WRAP_TO_REPLACE = 2;
            % replace potentially problematic indexes with Nan
            for i = 1: length ( BadIndexes )
                index_start = BadIndexes(i) -  NUM_SAMPLES_FROM_WRAP_TO_REPLACE ; 
                index_end = BadIndexes(i) +  NUM_SAMPLES_FROM_WRAP_TO_REPLACE ; 
                try
                    fixedPos( index_start : index_end ) = NaN;
                catch
                    fixedPos( BadIndexes(i) ) = NaN;
                end
            end


            % find start & end idxs of continuous chunks inbetween fictrac breaks
            breakIdxs = isnan(fixedPos);
            [chunks] = findContinuousData(breakIdxs);

            % correct for partial wrap between chunks --> "unwrap"
            for ch = 1:size(chunks,1)
                if ch ~= 1
                    lastChunkVal = fixedPos(chunks(ch - 1, 2));
                    while isnan(lastChunkVal)
                        lastChunkVal = fixedPos(chunks(ch - 1, 2) - 1);
                    end
                    firstChunkVal = fixedPos(chunks(ch, 1));
                    while isnan(firstChunkVal)
                        firstChunkVal = fixedPos(chunks(ch, 1) - 1);
                    end
                    chunkDiff = lastChunkVal - firstChunkVal; 
                    fixedPos(chunks(ch,1):chunks(ch,2)) = fixedPos(chunks(ch,1):chunks(ch,2)) + chunkDiff;
                end
            end

            unwrappedPos = fixedPos;
        catch
            disp("Couldn't fix")
        end

    end
end

% find indexes where the unwrapping happened (tolerance = pi)
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
    try
        cleanedPos(nanIDX) = cleanedPos(nanIDX - 1);
    catch
        cleanedPos(nanIDX) = cleanedPos(nanIDX + 1);
    end
    
    % find any remaining NaN
    nanIDX  = find( isnan(cleanedPos) );
end

% downsample position data to match FicTrac's output
[downsampled_cleanedPos,~] = resample_new(cleanedPos,(FicTracRate/2),sampleRate,[]); 

% smooth the position array
filteredPosition = smoothdata(downsampled_cleanedPos,'loess',10);

% convert to proper units
 if yaw == 1
        filteredPosition_unitConv = ( filteredPosition / (2*pi) ) * 360; % degrees
    else
        filteredPosition_unitConv = ( filteredPosition / (2*pi) ) * (pi*9); % mm
end


% differentiate into velocity
accumulatedPositionOut = resample_new(filteredPosition_unitConv,sampleRate,(FicTracRate/2),[]); 
velocity = gradient( filteredPosition ) .* (FicTracRate/2) ; 

 if yaw == 1
        velocity_unitConv = ( velocity / (2*pi) ) * 360; % degrees/sec
    else
        velocity_unitConv = ( velocity / (2*pi) ) * (pi*9); % mm/sec
end
 

% Calculate the distribution and take away values that are below 0.5% and above 99.5%
percentile25AV = prctile(velocity_unitConv,0.5);
percentile975AV = prctile(velocity_unitConv,99.5);
boundedVelocity = velocity_unitConv;
boundedVelocity(velocity_unitConv<percentile25AV | velocity_unitConv>percentile975AV) = NaN;

% Linearly interpolate to replace the NaNs with values.
[pointsVectorAV] = find(~isnan(boundedVelocity));
valuesVectorAV = boundedVelocity(pointsVectorAV);
xiAV = 1:length(boundedVelocity);
interpVel = interp1(pointsVectorAV,valuesVectorAV,xiAV);

if isnan(interpVel(1)) && ~isnan(interpVel(2))
    isnan(interpVel(1)) == isnan(interpVel(2)); %% urm, do better
end

% smooth velocity
velocityOut = resample_new(interpVel,sampleRate,(FicTracRate/2));
if yaw
    if ~isempty(BadIndexes) || setTo0
        try
            figure(8);clf; 
            subplot(4,1,1)
            plot(unwrap( posRadians ))
            hold on
            if isempty(startMovIdx)
                startMovIdx = length(posRadians); 
            end
            plot(startMovIdx,min((unwrap( posRadians ))),'ro')
            subplot(4,1,2)
            plot(abs(gradient(unwrap( posRadians ))))
            subplot(4,1,3)
            plot(fixedPos)
            subplot(4,1,4)
            plot(velocityOut)
        catch
            disp("couldn't display figure")
        end
    end
end
%% plotting to check how well unwrapping, cleaning and filtering worked
% can be commented out once you are happy with the parameters
 
% figure('Position',[50 50 800 900]),
% h(1) = subplot(8,1,1);
% plot(ficTracBallPosition)
% title('Raw voltage signal');
% ylim([0 10]);
% xlim([0 length(ficTracBallPosition)]);
% set(gca,'xticklabel',{[]})
% 
% h(2) = subplot(8,1,2);
% plot(posRadians)
% title('Signal in radians');
% ylim([0 2*pi]);
% xlim([0 length(posRadians)]);
% set(gca,'xticklabel',{[]})
% 
% h(3) = subplot(8,1,3);
% plot(cleanedPos)
% title('Unwrapped signal in radians');
% xlim([0 length(cleanedPos)]);
% set(gca,'xticklabel',{[]})
% 
% h(4) = subplot(8,1,4);
% plot(downsampled_cleanedPos)
% title('Downsampled signal in radians');
% xlim([0 length(downsampled_cleanedPos)]);
% set(gca,'xticklabel',{[]})
% 
% h(5) = subplot(8,1,5);
% plot(filteredPosition)
% title('lowpass Filtered position signal');
% xlim([0 length(filteredPosition)]);
% set(gca,'xticklabel',{[]})
% 
% h(6) = subplot(8,1,6);
% plot(filteredPosition_unitConv)
% title('lowpass Filtered position signal in mm or deg');
% xlim([0 length(filteredPosition_unitConv)]);
% set(gca,'xticklabel',{[]})
% 
% h(7) = subplot(8,1,7);
% plot(velocity)
% title('velocity signal in mm or deg');
% xlim([0 length(velocity)]);
% set(gca,'xticklabel',{[]})
% 
% h(8) = subplot(8,1,8);
% plot(velocityOut)
% title('Smoothed velocity signal in mm or deg');
% xlim([0 length(velocityOut)]);
% xlabel('Time (frames)');
% 
% linkaxes(h,'x');

end

