folders = get_folders_ephys(rootDir);

for f = 1:size(folders,1)
    folder = folders(f).folder;
    
    if strcmp(folder(end),'.')
                folder = folder(1:end-2); 
    end

    load(fullfile(folder,'pro_behaviourData.mat'))
    load(fullfile(folder,'pro_trialData.mat'))
    noPulse_behaviourData = processed_behaviourData{1}; 
    noPulse_trialData = processed_trialData{1}; 

    % find pulse indicies
    % plan: remove 10 seconds following 300 & 500ms pulses, 5 seconds following
    % 100 & 200 ms pulses
    try
        window = 10; 
        [pulse_table] = detect_iontoPulses(noPulse_behaviourData, window, 1000); 
        remove_idx = []; 

    for p = 1:size(pulse_table,1)
       pulse = pulse_table.pulseLength(p); 
       if pulse == 300 || pulse == 500
           pulse_idx = pulse_table.pulseStart(p):pulse_table.windowEnd(p); 
       else
           pulse_idx = pulse_table.pulseStart(p):pulse_table.windowEnd(p) - ((window * 1000)/2); 
       end

       remove_idx = [remove_idx, pulse_idx];
    end
    noPulse_behaviourData(remove_idx,:) = [];
    noPulse_trialData(remove_idx,:) = [];
    catch
        disp([folder, ': not a pulse trial'])
    end

    [trialGoal, ~, ~, ~] = plotHeading(noPulse_behaviourData, 60, 1000, folder,1);

    folders_cell{1} = folder;
    activityVSbehaviour_ionto_mean(folders_cell,noPulse_behaviourData,noPulse_trialData, 'fRate',1) 
    activityVSbehaviour_ionto_mean(folders_cell,noPulse_behaviourData,noPulse_trialData, 'Vm',1)

    prefHead = input('pref heading? ');

    try
        Ionto_behaviourVSactivity_lineplots(folder, prefHead,120,noPulse_trialData, noPulse_behaviourData,1)
    catch
        disp([folder, ': not enough headings sampled'])
    end
end