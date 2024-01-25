function [ velocityOut, accumulatedPosition ] = CalculateVelocity(accumulatedPosition, maxFlyVelocity, yaw, dsamp)
    ephysSettings
    
    accumulatedPosition = smoothdata(accumulatedPosition,'lowess',(100/1000)*settings.sampRate);
    
    if yaw == 1
        accumulatedPosition = ( accumulatedPosition / (2*pi) ) * 360;
    else
        accumulatedPosition = ( accumulatedPosition / (2*pi) ) * (pi*9);
    end
    
    if dsamp == 1
        accumulatedPosition_down = resample(accumulatedPosition,25,settings.sampRate);
    end
    velocityOut = gradient( accumulatedPosition_down ) .* 25 ;

    %% remove velocity values that are too large to be possible for the fly
   
    velocityOut(velocityOut > maxFlyVelocity) = maxFlyVelocity;
    velocityOut(velocityOut < -maxFlyVelocity) = -maxFlyVelocity;
    
        

end