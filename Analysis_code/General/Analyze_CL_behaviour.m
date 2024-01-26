function [disp_for, disp_side, disp_yaw, frX, frY, angle, vel_for, vel_side, vel_yaw] = Analyze_CL_behaviour(trial_data, pattern_num)

ephysSettings;

% When assigning fictrac directions, x will be forward, y will be side
ft_for = trial_data( :, settings.raw.x );
ft_yaw = trial_data( :, settings.raw.yaw);
ft_side = trial_data( :, settings.raw.y);
panelsX = trial_data( :, settings.raw.panelsX );
panelsY = trial_data( :, settings.raw.panelsY );
panels = [panelsX,panelsY];

%zero_vel_data = load( [experiment_dir '\' settings.zero_params_filename ] );

[ vel_for, disp_for ] = calcVel_analysis(ft_for, settings.sampRate, settings.fictracRate, 0); 
[ vel_yaw, disp_yaw ] = calcVel_analysis(ft_yaw, settings.sampRate, settings.fictracRate, 1);
[ vel_side, disp_side ] = calcVel_analysis(ft_side, settings.sampRate, settings.fictracRate, 0);

if pattern_num == 1 || pattern_num == 6
    yframeStep = 8;
elseif pattern_num == 3
    yframeStep = 24;
elseif pattern_num == 5
   yframeStep = 1;
end

[ frX, frY, angle] = process_panel_360( panels,settings.panels.numFramesX,settings.panels.numFramesY,yframeStep, pattern_num);

end

