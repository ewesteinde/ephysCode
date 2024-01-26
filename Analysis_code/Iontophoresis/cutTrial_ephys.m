folder = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\ionto_exps\meno_candidates\20220727_fly2_cell1\trial_2';
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

goal1_idx = find((noPulse_behaviourData.time > 210 & noPulse_behaviourData.time < 447)| (noPulse_behaviourData.time > 448 & noPulse_behaviourData.time < 600));

goal1_behaviourData =  noPulse_behaviourData(goal1_idx,:);
goal1_trialData =  noPulse_trialData(goal1_idx,:);
                    
    plotHeading(goal1_behaviourData, 60, 1000, folder,0);

        folders_cell{1} = folder;
        activityVSbehaviour_ionto_mean(folders_cell,goal1_behaviourData,goal1_trialData, 'fRate',0) 
        activityVSbehaviour_ionto_mean(folders_cell,goal1_behaviourData,goal1_trialData, 'Vm',0)

        prefHead = 180;

        try
            Ionto_behaviourVSactivity_lineplots(folder, prefHead,120,goal1_trialData, goal1_behaviourData,0)
        catch
            disp([folder, ': not enough headings sampled'])
        end
        
goal2_idx = find( (noPulse_behaviourData.time > 67 & noPulse_behaviourData.time < 221) | (noPulse_behaviourData.time > 448 & noPulse_behaviourData.time < 600));

goal2_behaviourData =  noPulse_behaviourData(goal2_idx,:);
goal2_trialData =  noPulse_trialData(goal2_idx,:);
                    
    [trialGoal, ~, ~, ~] = plotHeading(goal2_behaviourData, 60, 1000, folder,1);

        folders_cell{1} = folder;
        activityVSbehaviour_ionto_mean(folders_cell,goal2_behaviourData,goal2_trialData, 'fRate',0) 
        activityVSbehaviour_ionto_mean(folders_cell,goal2_behaviourData,goal2_trialData, 'Vm',0)

        prefHead = input('pref heading? ');

        try
            Ionto_behaviourVSactivity_lineplots(folder, prefHead,120,goal2_trialData, goal2_behaviourData)
        catch
            disp([folder, ': not enough headings sampled'])
        end

        