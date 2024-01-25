function closedLoop_yFunc(pattern, func, start_position)
%% begins closedLoop setting in panels
    freq = 50;
    Panel_com('stop');
    %set pattern number
    Panel_com('set_pattern_id', pattern);
    %set open loop for y & closed loop for x
    Panel_com('set_mode', [3, 4]); 
    Panel_com('set_funcy_freq' , freq);
    Panel_com('set_posFunc_id', [2, func]);
    Panel_com('set_position', [start_position, 1]);

    %quiet mode on
    Panel_com('quiet_mode_on');
end