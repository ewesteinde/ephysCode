function closeFicTrac()
    closed = 0; 
    while closed == 0
        try 
            [~, cmdOut] = system('tasklist | findstr /i "fictrac.exe"');
            cmdOut = strsplit(cmdOut);
            pid = cmdOut{2};
            system(['"C:\Code\windows-kill_x64_1.1.4_lib_release\windows-kill.exe" -SIGINT ', pid])
        catch
            disp('fictrac OL_CL is closed')
            closed = 1; 
        end 
    end
end