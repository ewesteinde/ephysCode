function SetTrialAveBehaviourFig(rootDir, pulse_window, sampRate)
    arguments
        rootDir string
        pulse_window double = 10
        sampRate double = 1000
    end
    
%% find appropriate folders
folders_all = get_folders_ephys_behaviour(rootDir, 1); 
count = 1; 
for f = 1:length(folders_all)
    if ~isempty(regexp(folders_all(f).folder, 'ionto_sCL')) || ~isempty(regexp(folders_all(f).folder, 'ionto_CL')) || ~isempty(regexp(folders_all(f).folder, 'ionto_jCL'))
        folders(count) = folders_all(f);
        count = count + 1;
    end
end


%% Process each folder
folderNum = length(folders);
fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
countFail = 1;
for ff = 1:folderNum
      %% Get folder information
      folder = folders(ff).folder;
      
    load(fullfile(folder,'pro_trialData.mat'));
    load(fullfile(folder,'pro_behaviourData.mat'));

    for t = 1:length(processed_behaviourData)
        close all
        [pulse_table] = detect_iontoPulses(processed_behaviourData{t}, pulse_window, sampRate);
        pulseLengths = unique(pulse_table.pulseLength);
        %% plot ave of equal sized pulses
        for p = 1:length(pulseLengths)
            pulse = pulseLengths(p);
            pulses = pulse_table(pulse_table.pulseLength == pulse,:);
            figure; clf;
            numPulses = size(pulses,1);
            
            
            time = [-pulse_window:1/sampRate:pulse_window];
            h(6) = subplot(16,1,1:3);
            try
                ScOAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                for i = 1:numPulses
                    ScOAve1 = ScOAve1 + processed_trialData{t}.scaledOutput(pulses.windowStart(i):pulses.windowEnd(i));
                    plot(time, processed_trialData{t}.scaledOutput(pulses.windowStart(i):pulses.windowEnd(i)))
                    hold on
                end
            catch
                ScOAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                for i = 1:numPulses -1
                    ScOAve1 = ScOAve1 + processed_trialData{t}.scaledOutput(pulses.windowStart(i):pulses.windowEnd(i));
                    plot(time, processed_trialData{t}.scaledOutput(pulses.windowStart(i):pulses.windowEnd(i)))
                    hold on
                end
            end
            ScOAve1 = ScOAve1/i;
            plot(time, ScOAve1, 'k', 'LineWidth',1.5) 
            ylabel('V (mV)')
            xlim([time(1),time(end)])


            %% angle 
            time = [-pulse_window:1/sampRate:pulse_window];
            h(1) = subplot(16,1,4:6);
            try
                for i = 1:numPulses
                    plot(time, processed_behaviourData{t}.angle(pulses.windowStart(i):pulses.windowEnd(i)))
                    hold on
                end
            catch
                for i = 1:numPulses -1
                    plot(time, processed_behaviourData{t}.angle(pulses.windowStart(i):pulses.windowEnd(i)))
                    hold on
                end
            end
            ylabel('cue angle')
            xlim([time(1),time(end)])

    %% for
            time = [-pulse_window:1/sampRate:pulse_window];
            h(2) = subplot(16,1,7:9);
            try
                forAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                for i = 1:numPulses
                    forAve1 = forAve1 + processed_behaviourData{t}.vel_for(pulses.windowStart(i):pulses.windowEnd(i));
                    plot(time, processed_behaviourData{t}.vel_for(pulses.windowStart(i):pulses.windowEnd(i)))
                    hold on
                end
            catch
                forAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                for i = 1:numPulses - 1
                    forAve1 = forAve1 + processed_behaviourData{t}.vel_for(pulses.windowStart(i):pulses.windowEnd(i));
                    plot(time, processed_behaviourData{t}.vel_for(pulses.windowStart(i):pulses.windowEnd(i)))
                    hold on
                end
            end
            forAve1 = forAve1/i;
            plot(time, forAve1, 'k', 'LineWidth',1.5) 
            ylabel('Vf (mm/sec)')
            xlim([time(1),time(end)])

    %% yaw
            time = [-pulse_window:1/sampRate:pulse_window];
            h(3) = subplot(16,1,10:12);
            try
                yawAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                for i = 1:numPulses
                    yawAve1 = yawAve1 + abs(processed_behaviourData{t}.vel_yaw(pulses.windowStart(i):pulses.windowEnd(i)));
                    plot(time, abs(processed_behaviourData{t}.vel_yaw(pulses.windowStart(i):pulses.windowEnd(i))))
                    hold on
                end
            catch
                yawAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                for i = 1:numPulses - 1
                    yawAve1 = yawAve1 + abs(processed_behaviourData{t}.vel_yaw(pulses.windowStart(i):pulses.windowEnd(i)));
                    plot(time, abs(processed_behaviourData{t}.vel_yaw(pulses.windowStart(i):pulses.windowEnd(i))))
                    hold on
                end
            end
            yawAve1 = yawAve1/i;
            plot(time, yawAve1, 'k', 'LineWidth',1.5) 
            ylabel('abs Vy (deg/sec)')
            xlim([time(1),time(end)])

    %% side

            time = [-pulse_window:1/sampRate:pulse_window];
            h(4) = subplot(16,1,13:15);
            try
                sideAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                for i = 1:numPulses
                    sideAve1 = sideAve1 + abs(processed_behaviourData{t}.vel_side(pulses.windowStart(i):pulses.windowEnd(i)));
                    plot(time, abs(processed_behaviourData{t}.vel_side(pulses.windowStart(i):pulses.windowEnd(i))))
                    hold on
                end
            catch
                sideAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                for i = 1:numPulses - 1
                    sideAve1 = sideAve1 + abs(processed_behaviourData{t}.vel_side(pulses.windowStart(i):pulses.windowEnd(i)));
                    plot(time, abs(processed_behaviourData{t}.vel_side(pulses.windowStart(i):pulses.windowEnd(i))))
                    hold on
                end
            end
            sideAve1 = sideAve1/i;
            plot(time, sideAve1, 'k', 'LineWidth',1.5) 
            ylabel('abs Vs (mm/sec)')
            xlim([time(1),time(end)])


            h(5) = subplot(16,1,16);
            plot(time, processed_behaviourData{t}.stim(pulses.windowStart(1):pulses.windowEnd(1)),'k')
            ylabel('Ionto Pulse')
            xlabel('Time (s)')
            sgtitle(strcat(num2str(pulse),'ms Pulse ',num2str(1),'-',num2str(numPulses), ' 10mM ATP'))
            xlim([time(1), time(end)])
            linkaxes(h,'x')

            saveas(gcf,fullfile(folder,'figures',[num2str(pulse),'msPulse_expAve.fig']))
        end
    end
end
end
