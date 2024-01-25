function down_tData = downSampleTDAQdata(trialData, fHz)

ephysSettings
%[y, t] = resample_new(x, fs_new, fs_old)

[current,time] = resample_new(trialData.current,fHz,(settings.sampRate));
[voltage,~] = resample_new(trialData.voltage,fHz,(settings.sampRate));
[scaledOutput,~] = resample_new(trialData.scaledOutput,fHz,(settings.sampRate));


down_tData = table(time, current, voltage, scaledOutput);

end