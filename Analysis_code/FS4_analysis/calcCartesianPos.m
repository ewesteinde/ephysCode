function [xPos, yPos] = calcCartesianPos(bData)

            %% Helen's path code

            %for nTrial = numTrials
                %notNan_idx = find(~isnan(bData.angle) & ~isnan(bData.vel_for) & ~isnan(bData.vel_side));
                yawAngPos = bData.angle;
                fwdAngVel = bData.vel_for;
                slideAngVel = bData.vel_side;

                sampRate = 1000; 
                % conversion factor between degrees and mm
                circum = 9 * pi; % circumference of ball, in mm
                mmPerDeg = circum / 360; % mm per degree of ball

                % position incorporating heading - as if fly were walking on x-y plane,
                %  x-y coordinates at each time point
                % start with fly at (0,0) and facing 0 deg
                zeroedYawAngPos = yawAngPos - yawAngPos(1); 

                % movement in x (in degrees) at each time point
                xChangePos = (fwdAngVel ./ sampRate) .* sind(zeroedYawAngPos) + ...
                    (slideAngVel ./ sampRate) .* sind(zeroedYawAngPos + 90);  

                % x position in mm (i.e. x-coordinate of fly's position at each time 
                %  point), starts at 0
                xPos = (cumsum(xChangePos) - xChangePos(1)) .* mmPerDeg;


                % movement in y (in degrees) at each time point
                yChangePos = (fwdAngVel ./ sampRate) .* cosd(zeroedYawAngPos) + ...
                    (slideAngVel ./ sampRate) .* cosd(zeroedYawAngPos + 90);

                % y position in mm (i.e. y-coordinate of fly's position at each time 
                %  point), starts at 0
                yPos = (cumsum(yChangePos) - yChangePos(1)) .* mmPerDeg;

                
end
%close all 