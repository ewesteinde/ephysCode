function configFictrac()

expDir = 'C:\Code\FicTrac\sample';
FT_PATH = 'C:\Code\FicTrac\bin\Release\configGui.exe';
cmdStr = ['cd "', expDir, '" & start "" "',  FT_PATH, ...
        '" config_current.txt & exit'];
system(cmdStr);
end