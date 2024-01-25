function runClosedLoopPattern()
ephysSettings
Panel_com('stop');
Panel_com('all_off');

runFictrac()

% run CL pattern in between recorded trials
pattern = settings.panels.barPattern;
startPosition = 1; 
closedLoop(pattern, startPosition)
Panel_com('start');
end