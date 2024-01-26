function runPanelsClosedLoop()
    % close fictrac if it's already running & stop panels
    closed = 0; 
    while closed == 0
        try 
            [~, cmdOut] = system('tasklist | findstr /i "fictrac.exe"');
            cmdOut = strsplit(cmdOut);
            pid = cmdOut{2};
            system(['"C:\Code\windows-kill_x64_1.1.4_lib_release\windows-kill.exe" -SIGINT ', pid])
        catch
            disp('fictrac CL is closed')
            closed = 1; 
        end 
    end

    Panel_com('stop')

    % Start FicTrac in background from current experiment directory (config file must be in directory)
    expDir = 'C:\Code\FicTrac\sample';
    FT_PATH = 'C:\Code\FicTrac\bin\Release\fictrac.exe';
    cmdStr = ['cd "', expDir, '" & start "" "',  FT_PATH, ...
            '" config_current.txt & exit'];
    system(cmdStr);

    ephysSettings
    % run CL pattern in between recorded trials
    pattern = settings.panels.barPattern;
    startPosition = 1; 
    closedLoop(pattern, startPosition)
    Panel_com('start');

    closedLoop(settings.panels.barPattern, 1)
    Panel_com('start');
    
 % Call socket_client_360 to open socket connection from fictrac to Phiget22 device
    Socket_PATH = 'C:\Code\FicTrac\scripts';
    SOCKET_SCRIPT_NAME = 'socket_client_360.py';
    cmdstring = ['cd "' Socket_PATH '" & python ' SOCKET_SCRIPT_NAME ' &'];
    [status] = system(cmdstring, '-echo');
end