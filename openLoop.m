function openLoop(pattern, func, start_position)
%% begins closedLoop setting in panels
    freq = 50;
    Panel_com('stop');
    %set pattern number
    Panel_com('set_pattern_id', pattern);
    %set open loop for x
    Panel_com('set_mode', [4, 0]);
    Panel_com('set_funcX_freq' , freq);
    Panel_com('set_posFunc_id', [1, func]);
    Panel_com('set_position', [start_position, 1]);
    %quiet mode on
    Panel_com('quiet_mode_on');
end