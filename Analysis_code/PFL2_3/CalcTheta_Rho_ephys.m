function [rho, theta] = CalcTheta_Rho_ephys(behaviourData,minVel,window,data_idx)


if strcmp(data_idx,'all')
    data_idx = 1:length(behaviourData.angle);
end

angle_temp = resample_new(behaviourData.angle(data_idx),100,1000,1);
angle_temp(angle_temp > 180) = 180; 
angle_temp(angle_temp < -180) = -180; 
vf = resample_new(behaviourData.vel_for(data_idx),100,1000,0);
vy = resample_new(behaviourData.vel_yaw(data_idx),100,1000,0);
vs = resample_new(behaviourData.vel_side(data_idx),100,1000,0);
count = 1;
window = window * 100; 
step = 1;
for i = 1 - window/2:step:length(angle_temp) - window/2
    idx = i:i + window; 

    if idx(end) > length(angle_temp) && idx(1) < 1 
        idx = 1:1:length(angle_temp);
    elseif idx(end) > length(angle_temp)
        idx = idx(1):1:length(angle_temp);
    elseif idx(1) < 1 
        idx = 1:idx(end); 
    end

    speed = abs(vf) + abs(vs) + abs(deg2rad(vy)*4.5);
    angles_flyFor = angle_temp(speed > minVel); 
        angles_flyFor = angle_temp(idx);
        if ~isempty(angles_flyFor)
            x = cosd(angles_flyFor); 
            y = sind(angles_flyFor); %my arena has - angles to the left of the fly, + to the right, multiply y component by -1 to align physical arena coords to polar plot angles
            idx_windows{count,1} = idx;
            mean_headingVectors(1,count)= sum(x)/length(x); 
            mean_headingVectors(2,count)= sum(y)/length(y);
            count = count + 1;
        else
            mean_headingVectors(1,count)= nan; 
            mean_headingVectors(2,count)= nan; 
            idx_windows{count,1} = idx;
            count = count + 1;
        end
end

 rho = sqrt(mean_headingVectors(1,:).^2 + mean_headingVectors(2,:).^2);
 rho = resample_new(rho,1000,100,0);
 theta = atan2(mean_headingVectors(2,:),mean_headingVectors(1,:));
 theta = resample_new(theta,1000,100,0); 