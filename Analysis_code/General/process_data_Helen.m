function fictrac = process_data_Helen(trial_data, pattern_num)

ephysSettings;

% When assigning fictrac directions, x will be forward, y will be side
rawDAQ.for = trial_data( :, settings.raw.x );
rawDAQ.yaw = trial_data( :, settings.raw.yaw);
rawDAQ.side = trial_data( :, settings.raw.y);
numSamp = size(rawDAQ.for,1);
daqTime = (0:(numSamp-1))/settings.sampRate;
BALL_DIAM = 9; 

panelsX = trial_data( :, settings.raw.panelsX );
panelsY = trial_data( :, settings.raw.panelsY );
panels = [panelsX,panelsY];

% conversion factor between degrees and mm
circum = BALL_DIAM * pi; % circumference of ball, in mm
mmPerDeg = circum / 360; % mm per degree of ball

[ vel_for_ang, disp_for_ang ] = ficTracSignalDecoding_helen(rawDAQ.for, settings.sampRate, settings.fictracRate/2,2500); 
     
vel_for = vel_for_ang .* mmPerDeg; % velocity in mm/sec

% cumulative forward position in mm, where start of trial is at 0
disp_for = (disp_for_ang - disp_for_ang(1)) .* mmPerDeg; 
    
[ vel_yaw, disp_yaw_ang ] = ficTracSignalDecoding_helen(rawDAQ.yaw, settings.sampRate, settings.fictracRate/2,2500);

% wrap yaw angular position to 360 deg instead of it being cumulative
disp_yaw = wrapTo360(disp_yaw_ang);

[ vel_side_ang, disp_side_ang ] = ficTracSignalDecoding_helen(rawDAQ.side, settings.sampRate, settings.fictracRate/2,2500);

 vel_side = vel_side_ang .* mmPerDeg; % velocity in mm/sec
% cumulative slide position in mm, where start of trial is at 0
disp_side = (disp_side_ang-disp_side_ang(1)) .* mmPerDeg;    


 % position incorporating heading - as if fly were walking on x-y plane,
    %  x-y coordinates at each time point
    % start with fly at (0,0) and facing 0 deg
    zeroedYawAngPos = disp_yaw_ang - disp_yaw_ang(1); 
    
    % movement in x (in degrees) at each time point
    xChangePos = (vel_for_ang ./ settings.sampRate) .* sind(zeroedYawAngPos) + ...
        (vel_side_ang ./ settings.sampRate) .* sind(zeroedYawAngPos + 90);  

    % x position in mm (i.e. x-coordinate of fly's position at each time 
    %  point), starts at 0
    xPos = (cumsum(xChangePos) - xChangePos(1)) .* mmPerDeg;
   
    % movement in y (in degrees) at each time point
    yChangePos = (vel_for_ang ./ settings.sampRate) .* cosd(zeroedYawAngPos) + ...
        (vel_side_ang ./ settings.sampRate) .* cosd(zeroedYawAngPos + 90);

    % y position in mm (i.e. y-coordinate of fly's position at each time 
    %  point), starts at 0
    yPos = (cumsum(yChangePos) - yChangePos(1)) .* mmPerDeg;



[fictrac.vel_for,time] = resample_new(vel_for, 1000, settings.sampRate);
[fictrac.disp_for,~] = resample_new(disp_for, 1000, settings.sampRate);
[fictrac.vel_yaw,~] = resample_new(vel_yaw, 1000, settings.sampRate);
[fictrac.disp_yaw,~] = resample_new(disp_yaw, 1000, settings.sampRate);
[fictrac.vel_side,~] = resample_new(vel_side, 1000, settings.sampRate);
[fictrac.disp_side,~] = resample_new(disp_side, 1000, settings.sampRate);
fictrac.xPos = xPos;
fictrac.yPos = yPos;
fictrac.daqTime = daqTime; 
fictrac.downTime = time; 



if pattern_num == 1
    yframeStep = 8;
elseif pattern_num == 3
    yframeStep = 24;
elseif pattern_num == 5
   yframeStep = 1;
end
[ frX, frY, angle] = process_panel_360( panels,settings.panels.numFramesX,settings.panels.numFramesY,yframeStep, pattern_num);

fictrac.frX = frX;
fictrac.frY = frY;
fictrac.angle = angle; 


end