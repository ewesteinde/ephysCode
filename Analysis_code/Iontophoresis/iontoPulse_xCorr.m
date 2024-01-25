function iontoPulse_xCorr(rootDir, pulse_window, sampRate)
    arguments
        rootDir string
        pulse_window double = 10
        sampRate double = 1000
    end
    
%% find appropriate folders
folders_all = get_folders_ephys_behaviour(rootDir, 1); 
count = 1; 

if isempty(regexp(rootDir,'topTrials'))
    for f = 1:length(folders_all)
        if ~isempty(regexp(folders_all(f).folder, 'ionto_sCL')) %|| ~isempty(regexp(folders_all(f).folder, 'ionto_CL')) || ~isempty(regexp(folders_all(f).folder, 'ionto_jCL'))
            folders(count) = folders_all(f);
            count = count + 1;
        end
    end
else
    folders = folders_all;
end


%% Process each folder
folderNum = length(folders);
fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
countFail = 1;
pulse_table = [];
for ff = 1:folderNum
      %% Get folder information
      folder = folders(ff).folder;
      
    load(fullfile(folder,'pro_trialData.mat'));
    all_trialData{ff} = processed_trialData{1};
    
    load(fullfile(folder,'pro_behaviourData.mat'));
    all_behaviourData{ff} = processed_behaviourData{1};
    
    [trialPulses] = detect_iontoPulses(all_behaviourData{ff}, pulse_window, sampRate);
    pulseLengths = unique(trialPulses.pulseLength);
    pulse_table{ff} = trialPulses;
end

        %% plot ave of equal sized pulses
    grey = [0.75,0.75,0.75];
    red = [0.75, 0, 0, 0.25];
    corrLag = 4; 
    
    for p = 1:length(pulseLengths)
        forCorr_Ave = zeros( corrLag*sampRate*2 + 1,1);
        yawCorr_Ave = zeros( corrLag*sampRate*2 + 1,1);
        sideCorr_Ave = zeros( corrLag*sampRate*2 + 1,1);
        fig = figure; clf;
        h1 = subplot(3,1,1);
        h2 = subplot(3,1,2);
        h3 = subplot(3,1,3);
        hold([h1,h2,h3],'on')
        for ff = 1:folderNum
        pulse = pulseLengths(p);
        pulses = pulse_table{ff}(pulse_table{ff}.pulseLength == pulse,:);
        numPulses = size(pulses,1);
        
        if numPulses > 0

            try
                vfCorr_trialAve = zeros(corrLag*sampRate*2 + 1,1);
                for i = 1:numPulses
                    [vfCorr, vfLags] = xcorr(all_trialData{ff}.scaledOutput(pulses.windowStart(i):pulses.windowEnd(i)),all_behaviourData{ff}.vel_for(pulses.windowStart(i):pulses.windowEnd(i)),corrLag*sampRate,'normalized');
                    plot(h1,vfLags*1/sampRate, vfCorr,'color',grey)
                    vfCorr_trialAve = vfCorr_trialAve + vfCorr; 
                end
            catch
                vfCorr_trialAve = zeros(corrLag*sampRate*2 + 1,1);
                for i = 1:numPulses -1
                    [vfCorr, vfLags] = xcorr(all_trialData{ff}.scaledOutput(pulses.windowStart(i):pulses.windowEnd(i)),all_behaviourData{ff}.vel_for(pulses.windowStart(i):pulses.windowEnd(i)),corrLag*sampRate,'normalized');
                    plot(h1,vfLags*1/sampRate, vfCorr,'color',grey)
                    vfCorr_trialAve = vfCorr_trialAve + vfCorr; 
                end
            end
            vfCorr_trialAve = vfCorr_trialAve/i;
            plot(h1,vfLags*1/sampRate,vfCorr_trialAve,'k','LineWidth',1.5)


             %% angle 
%             time = [-pulse_window:1/sampRate:pulse_window];
%             try
%                 for i = 1:numPulses
%                     plot(h1,time, all_behaviourData{ff}.angle(pulses.windowStart(i):pulses.windowEnd(i)))
%                     hold on
%                 end
%             catch
%                 for i = 1:numPulses -1
%                     plot(h1,time, all_behaviourData{ff}.angle(pulses.windowStart(i):pulses.windowEnd(i)))
%                     hold on
%                 end
%             end
%             ylabel(h1,'cue angle')
%             xlim(h1,[time(1),time(end)])

    %% for
            time = [-pulse_window:1/sampRate:pulse_window];

            try
                forAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                for i = 1:numPulses
                    forAve1 = forAve1 + all_behaviourData{ff}.vel_for(pulses.windowStart(i):pulses.windowEnd(i));
                    plot(h2,time, all_behaviourData{ff}.vel_for(pulses.windowStart(i):pulses.windowEnd(i)),'color',grey)
                    hold on
                end
            catch
                forAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                for i = 1:numPulses - 1
                    forAve1 = forAve1 + all_behaviourData{ff}.vel_for(pulses.windowStart(i):pulses.windowEnd(i));
                    plot(h2,time, all_behaviourData{ff}.vel_for(pulses.windowStart(i):pulses.windowEnd(i)),'color',grey)
                    hold on
                end
            end
            forAve1 = forAve1/i;
            ylabel(h2,'Vf (mm/sec)')
            xlim(h2,[time(1),time(end)])

    %% yaw
            time = [-pulse_window:1/sampRate:pulse_window];

            try
                yawAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                for i = 1:numPulses
                    yawAve1 = yawAve1 + abs(all_behaviourData{ff}.vel_yaw(pulses.windowStart(i):pulses.windowEnd(i)));
                    plot(h3,time, abs(all_behaviourData{ff}.vel_yaw(pulses.windowStart(i):pulses.windowEnd(i))),'color',grey)
                    hold on
                end
            catch
                yawAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                for i = 1:numPulses - 1
                    yawAve1 = yawAve1 + abs(all_behaviourData{ff}.vel_yaw(pulses.windowStart(i):pulses.windowEnd(i)));
                    plot(h3,time, abs(all_behaviourData{ff}.vel_yaw(pulses.windowStart(i):pulses.windowEnd(i))),'color',grey)
                    hold on
                end
            end
            yawAve1 = yawAve1/i;
            ylabel(h3,'abs Vy (deg/sec)')
            xlim(h3,[time(1),time(end)])

    %% side

            time = [-pulse_window:1/sampRate:pulse_window];

            try
                sideAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                for i = 1:numPulses
                    sideAve1 = sideAve1 + abs(all_behaviourData{ff}.vel_side(pulses.windowStart(i):pulses.windowEnd(i)));
                    plot(h4,time, abs(all_behaviourData{ff}.vel_side(pulses.windowStart(i):pulses.windowEnd(i))),'color',grey)
                    hold on
                end
            catch
                sideAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                for i = 1:numPulses - 1
                    sideAve1 = sideAve1 + abs(all_behaviourData{ff}.vel_side(pulses.windowStart(i):pulses.windowEnd(i)));
                    plot(h4,time, abs(all_behaviourData{ff}.vel_side(pulses.windowStart(i):pulses.windowEnd(i))),'color',grey)
                    hold on
                end
            end
            sideAve1 = sideAve1/i;
            ylabel(h4,'abs Vs (mm/sec)')
            xlim(h4,[time(1),time(end)])
            xlabel(h4,'Time (s)')

            ScOAve = ScOAve + ScOAve1;
            forCorr_Ave = forCorr_Ave + forAve1;
            yawCorr_Ave = yawCorr_Ave + yawAve1;
            sideCorr_Ave = sideCorr_Ave + sideAve1;
        end
        end

            ScOAve = ScOAve/ff;
            forCorr_Ave = forCorr_Ave/ff;
            yawCorr_Ave = yawCorr_Ave/ff;
            sideCorr_Ave = sideCorr_Ave/ff;

            plot(h6,time, ScOAve, 'LineWidth',1.5,'Color','k')
            plot(h2,time, forCorr_Ave, 'LineWidth',1.5,'Color','k')
            plot(h3,time, yawCorr_Ave, 'LineWidth',1.5,'Color','k')
            plot(h4,time, sideCorr_Ave, 'LineWidth',1.5,'Color','k')
            
            ylim(h6,[-80 -30])
            ylim(h2,[-5 12])
            ylim(h3,[0 300])
            ylim(h4, [0 8])
            
            rectangle(h2,'Position', [0, min(ylim(h2)), pulse/1000, max(ylim(h2)) + abs(min(ylim(h2)))], 'FaceColor',red, 'EdgeColor', 'none');
            rectangle(h3,'Position', [0, min(ylim(h3)), pulse/1000, max(ylim(h3)) + abs(min(ylim(h3)))], 'FaceColor',red, 'EdgeColor', 'none');
            rectangle(h4,'Position', [0, min(ylim(h4)), pulse/1000, max(ylim(h4)) + abs(min(ylim(h4)))], 'FaceColor',red, 'EdgeColor', 'none');
            rectangle(h6,'Position', [0, min(ylim(h6)), pulse/1000, max(ylim(h6)) + abs(min(ylim(h6)))], 'FaceColor',red, 'EdgeColor', 'none');
            
            sgtitle(strcat(num2str(pulse),'ms Pulse'))
            linkaxes([h2,h3,h4,h6],'x')
            set(gcf,'color','w')
            hallbutlast = [h2,h3,h6];
            set(hallbutlast,'XTickLabel',[])
            %set(h6,'YTickLabel',[])
            set(h6, 'YColor','k')
            set(hallbutlast, 'XColor','none')
            box off
            set(fig,'units','normalized','position',[0.25 0.15 0.25 0.75])
            set(fig,'Renderer','painters')
            %subplotsqueezeV(fig, 1.1)
            saveas(fig,fullfile(rootDir,[num2str(pulse),'msPulse_expAve.svg']))
    end
end

