function closedLoop(pattern, startPosition)
    %% begins closedLoop setting in panels
    Panel_com('stop');
    %set arena
    pause(.05)
    Panel_com('set_config_id', 1);
    %set brightness level
    pause(.05)
%     Panel_com('g_level_1');
    %set pattern number
    Panel_com('set_pattern_id', pattern);
    Panel_com('set_position', [96, startPosition]);
    %% 
    %set closed loop for x
    pause(.05)
    Panel_com('set_mode', [3, 0]);
    pause(.05)
    Panel_com('quiet_mode_on');
end