function down_bData = downSampleBDAQdata(behaviorData, fHz)

ephysSettings
%[y, t] = resample_new(x, fs_new, fs_old)
if length(behaviorData.vel_for) == length(behaviorData.frY)
    [vel_for,~] = resample_new(behaviorData.vel_for,fHz,(settings.sampRate));
    [vel_yaw,~] = resample_new(behaviorData.vel_yaw,fHz,(settings.sampRate));
    [vel_side,~] = resample_new(behaviorData.vel_side,fHz,(settings.sampRate));
else
    vel_for = behaviorData.vel_for; 
    vel_yaw = behaviorData.vel_yaw;
    vel_side = behaviorData.vel_side;
end

if length(behaviorData.disp_for) == length(behaviorData.frY)
    [disp_for,~] = resample_new(behaviorData.disp_for,fHz,(settings.sampRate));
    [disp_yaw,~] = resample_new(behaviorData.disp_yaw,fHz,(settings.sampRate));
    [disp_side,~] = resample_new(behaviorData.disp_side,fHz,(settings.sampRate));
else
    disp_for = behaviorData.disp_for; 
    disp_yaw = behaviorData.disp_yaw;
    disp_side = behaviorData.disp_side;
end

[angle,~] = resample_new(behaviorData.angle,fHz,(settings.sampRate));
angle = smoothdata(angle,'movmedian',20);
angle(angle > 180) = 180;
angle(angle < -180) = -180; 
%     int_angle = int_angle';

[frY,time] = resample_new(behaviorData.frY,fHz,(settings.sampRate));
[~, jumps,frY_clean] = detect_jumps_ephys(frY, 10,10, fHz);
frY = frY_clean;

down_bData = table(time, disp_for, disp_side, disp_yaw, angle, vel_for, vel_side, vel_yaw, frY, jumps);

end