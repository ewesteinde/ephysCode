function [ velocityOut, accumulatedPosition ] = ConvertUnits(accumulatedPosition, velocity, yaw)
    ephysSettings
    
%Calculate the distribution and take away values that are below 2.5% and above 97.5%

percentile25AV = prctile(velocity,2.5);
percentile975AV = prctile(velocity,97.5);
boundedVelocity = velocity;
boundedVelocity(velocity<percentile25AV |velocity>percentile975AV) = NaN;

% 8)Linearly interpolate to replace the NaNs with values.
[pointsVectorAV] = find(~isnan(boundedVelocity));
valuesVectorAV = boundedVelocity(pointsVectorAV);
xiAV = 1:length(boundedVelocity);
interpVel = interp1(pointsVectorAV,valuesVectorAV,xiAV);

% 9)Smooth
angularVelocity = smoothdata(interpVel,'rlowess',15);
    
    if yaw == 1
        accumulatedPosition = ( accumulatedPosition / (2*pi) ) * 360;
        velocityOut = ( velocity / (2*pi) ) * 360;
    else
        accumulatedPosition = ( accumulatedPosition / (2*pi) ) * (pi*9);
        velocityOut = ( velocity / (2*pi) ) * (pi*9);
    end

    %% remove velocity values that are too large to be possible for the fly
%    
%     velocityOut(velocityOut > maxFlyVelocity) = maxFlyVelocity;
%     velocityOut(velocityOut < -maxFlyVelocity) = -maxFlyVelocity;
%     
        

end