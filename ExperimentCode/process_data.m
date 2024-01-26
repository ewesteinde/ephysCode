function [disp_for, disp_side, disp_yaw, frX, frY, angle, vel_for, vel_side, vel_yaw] = process_data(trial_data, pattern_num)

ephysSettings;

% When assigning fictrac directions, x will be forward, y will be side
ft_for = trial_data( :, settings.raw.x );
ft_yaw(:,1) = trial_data( :, settings.raw.yaw);
ft_side = trial_data( :, settings.raw.y);
panelsX = trial_data( :, settings.raw.panelsX );
panelsY = trial_data( :, settings.raw.panelsY );
%panelsYgain = trial_data( :, settings.raw.yawgain);
panels = [panelsX,panelsY];

%zero_vel_data = load( [experiment_dir '\' settings.zero_params_filename ] );

[ vel_for, disp_for ] = ficTracSignalDecoding(ft_for, settings.sampRate, settings.fictracRate, 0); 

try
    ft_yaw(:,2) = disp_for;
catch
    ft_yaw(:,2) = disp_for(1:length(ft_yaw(:,1)));
end

[ vel_yaw, disp_yaw ] = ficTracSignalDecoding(ft_yaw, settings.sampRate, settings.fictracRate, 1);
[ vel_side, disp_side ] = ficTracSignalDecoding(ft_side, settings.sampRate, settings.fictracRate, 0);

if pattern_num == 1
    yframeStep = 8;
    settings.panels.numFramesY = 16;
elseif pattern_num == 3 || pattern_num == 2
    yframeStep = 24;
    settings.panels.numFramesY = 4;
elseif pattern_num == 5
   yframeStep = 1;
   settings.panels.numFramesY = 96;
end

[ frX, frY, angle] = process_panel_360( panels,settings.panels.numFramesX,settings.panels.numFramesY,yframeStep, pattern_num);

end


