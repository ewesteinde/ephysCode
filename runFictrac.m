function runFictrac()

expDir = 'C:\Code\FicTrac\sample';
FT_PATH = 'C:\Code\FicTrac\bin\Release\fictrac.exe';
cmdStr = ['cd "', expDir, '" & start "" "',  FT_PATH, ...
        '" config_current.txt & exit'];
system(cmdStr);


 % Call socket_client_360 to open socket connection from fictrac to Phiget22 device
Socket_PATH = 'C:\Code\FicTrac\scripts';
SOCKET_SCRIPT_NAME = 'socket_client_360.py';
cmdstring = ['cd "' Socket_PATH '" & python ' SOCKET_SCRIPT_NAME ' &'];
[status] = system(cmdstring, '-echo');

end