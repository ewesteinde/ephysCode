function FindMovementTransitions(rootDir)

folders = get_folders_ephys(rootDir); 


    for f = 1:size(folders,1)
        folder = folders(f).folder; 
        trialIdx = regexp(folder,'_trial');
        if ~isempty(trialIdx)
            [~,cellIdx] = regexp(folder,'cell_');
            folder_new = folder(1:cellIdx + 1);
            flies(f) = string(folder_new);
        else
            flies(f) = string(folder); 
        end    
    end

flies = flies';
allFlies = unique(flies);

%%
transitionThreshold = round(1.5,1);
runThreshold = round(1,1); 
bufferTime = 2;
next_stopLength = bufferTime * 1000; 
prev_startLength = bufferTime * 1000; 


for t = 1:length(allFlies)
    
    fly = allFlies(t); 
    fly_folders = folders(flies == fly);
    
    acount_stop = 1; 
    pcount_stop = 1;
    acount_start = 1; 
    pcount_start = 1; 
    for f = 1:size(fly_folders,1)


            folder = fly_folders(f).folder; 
            if strcmp(folder(end),'.')
                folder = folder(1:end-2); 
            end

            processedDir = fullfile(folder,'processedData');
            load(fullfile(processedDir,'pro_behaviourData.mat'),'pro_behaviourData')
            load(fullfile(processedDir,'pro_trialData.mat'),'pro_trialData')

            bData = pro_behaviourData{1};
            tData = pro_trialData{1}; 
            total_mov_mm = round(abs(bData.vel_for + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw)*4.5)),1);
            total_mov_smooth = smoothdata(total_mov_mm,'gaussian',1000); 
            start_transitions = zeros(size(total_mov_mm)); 
            prev_startIdx = 1; 
            
%             figure()
%             plot(bData.time,total_mov_mm)
%             hold on
%             plot(bData.time,total_mov_smooth)


            for i = 2:length(total_mov_mm)-1
                prev_startWindow = (i - prev_startLength):i-1;
                if prev_startWindow(1) < 1
                    prev_startWindow = 1:prev_startWindow(end);
                end

                next_stopWindow = i+1:i+next_stopLength;
                if next_stopWindow(end) > length(total_mov_mm)
                    next_stopWindow = next_stopWindow(1):length(total_mov_mm); 
                end

                if  (total_mov_mm(i) == transitionThreshold) && (all(total_mov_mm(prev_startWindow) < transitionThreshold)) && (all(total_mov_smooth(next_stopWindow) > runThreshold)) %&& (all(total_mov_mm(i + prev_startLength) > runThreshold))
                    start_transitions(i) = 1;
                    prev_startIdx = i; 
                end
            end

            stop_transitions = zeros(size(total_mov_mm)); 
            next_stopIdx = 2; 


            for i = 2:length(total_mov_mm)-1
                next_stopWindow = i+1:i+next_stopLength;
                if next_stopWindow(end) > length(total_mov_mm)
                    next_stopWindow = next_stopWindow(1):length(total_mov_mm); 
                end

                prev_startWindow = (i - prev_startLength):i-1;
                if prev_startWindow(1) < 1
                    prev_startWindow = 1:prev_startWindow(end);
                end

                if  (total_mov_mm(i) == transitionThreshold) && (all(total_mov_mm(next_stopWindow) < transitionThreshold)) && (all(total_mov_smooth(prev_startWindow) > runThreshold)) 
                    stop_transitions(i) = 1;
                    prev_stopIdx = i; 
                end
            end
% 
%             figure();
%             plot(bData.time,total_mov_mm)
%             hold on
%             plot(bData.time(logical(stop_transitions)),total_mov_mm(logical(stop_transitions)),'ko')
%             plot(bData.time(logical(start_transitions)),total_mov_mm(logical(start_transitions)),'ro')

            startIdx = find(start_transitions == 1); 
            stopIdx = find(stop_transitions == 1);

            activity_stop = []; 
            activity_start = [];

            transitionWindow = 2 * 1000; 
            windowTime = -2:1/1000:2;
            count = 1; 
            for start = 1:length(startIdx)
                startWindow = startIdx(start) - transitionWindow: startIdx(start) + transitionWindow;
                if startWindow(end) < length(bData.angle) +1 && startWindow(1) > 0 
                    try
                        activity_start(count,:) = tData.smoothVm(startWindow);
                    catch
                        activity_start(count,:) = tData.smooth_Vm(startWindow);
                    end
                    count = count + 1; 
                end    
            end

            transitionWindow = 2 * 1000; 
            windowTime = -2:1/1000:2;
            count = 1; 
            for stop = 1:length(stopIdx)
                stopWindow = stopIdx(stop) - transitionWindow: stopIdx(stop) + transitionWindow;
                if stopWindow(end) < length(bData.angle) +1 && stopWindow(1) > 0 
                    try
                        activity_stop(count,:) = tData.smoothVm(stopWindow);
                    catch
                        activity_stop(count,:) = tData.smooth_Vm(stopWindow);
                    end
                    total_mov_stop(count,:) = total_mov_mm(stopWindow);
    %                 [corr, lags] = xcorr(activity_stop(count,:),total_mov_stop(count,:),1000);
    %                 
    %                 stop_corr(count,:) = corr; 

                    count = count + 1; 
                end       
            end
% 
%             figure(); plot(windowTime,mean(activity_stop,1),'k','LineWidth',1.5)
%             xline(0)
%             figure(); plot(windowTime,mean(activity_start,1),'k','LineWidth',1.5)

            prefHead = calcPrefHead(bData, tData);

            for i = 1:size(activity_start,1)
                angleStart = startIdx(i) - (bufferTime * 1000) + 1;
                if angleStart < 1
                    angleStart = 1; 
                end
                angleEnd = startIdx(i) - 1;
                startAngle = circ_mean(deg2rad(bData.angle(angleStart:angleEnd)));

                prefHeadDiff = abs(rad2deg(angdiff(startAngle, deg2rad(prefHead))));

                if prefHeadDiff <= 45
                    prefStarts{t}(pcount_start,:) = smoothdata(activity_start(i,:),'movmean',75);
                    pcount_start = pcount_start + 1;
                elseif prefHeadDiff >= 135
                    antiPrefStarts{t}(acount_start,:) = smoothdata(activity_start(i,:),'movmean',75); 
                    acount_start = acount_start + 1;  
                end
            end

            for i = 1:size(activity_stop,1)
                angleStart = stopIdx(i) + 1;
                angleEnd = stopIdx(i) + 1 + (bufferTime * 1000);
                if angleEnd > size(total_mov_mm,1)
                    angleEnd = size(total_mov_mm,1); 
                end
                 startAngle = circ_mean(deg2rad(bData.angle(angleStart:angleEnd)));

                prefHeadDiff = abs(rad2deg(angdiff(startAngle, deg2rad(prefHead))));

                if prefHeadDiff <= 45
                    prefStops{t}(pcount_stop,:) = smoothdata(activity_stop(i,:),'movmean',75);
                    pcount_stop = pcount_stop + 1;
                elseif prefHeadDiff >= 135
                    antiPrefStops{t}(acount_stop,:) = smoothdata(activity_stop(i,:),'movmean',75);
                    acount_stop = acount_stop + 1;  
                end
            end
    end
    
            breaks = regexp(fly,'\');    
            temp = char(fly);
            
            SEM_prefStops = std(prefStops{t},[],1,'omitnan') / sqrt(size(prefStops{t},1));
            SEM_antiPrefStops = std(antiPrefStops{t},[],1,'omitnan') / sqrt(size(antiPrefStops{t},1));
            SEM_prefStarts = std(prefStarts{t},[],1,'omitnan') / sqrt(size(prefStarts{t},1));
            SEM_antiPrefStarts = std(antiPrefStarts{t},[],1,'omitnan') / sqrt(size(antiPrefStarts{t},1));
    
            meanPrefStops = mean(prefStops{t},1);
            meanantiPrefStops = mean(antiPrefStops{t},1);
            meanPrefStarts = mean(prefStarts{t},1);
            meanantiPrefStarts = mean(antiPrefStarts{t},1);
            
%             figure();
%             set(gcf,'color','w')
%             set(gcf,'renderer','painters')
%             sgtitle(temp(breaks(end-2) + 1:end))
%             keepIndex = ~isnan(SEM_prefStops);
%             SEMhigh = [meanPrefStops(keepIndex) + SEM_prefStops(keepIndex)]; 
%             SEMlow = [meanPrefStops(keepIndex) - SEM_prefStops(keepIndex)];
%             ax1 = subplot(2,1,1);
%             plot(windowTime, meanPrefStops,'k','LineWidth',1.5)
%             patch([windowTime(keepIndex) fliplr(windowTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%             title(ax1,'prefStops')
%             xline(ax1,0)
%             box off
%             keepIndex = ~isnan(SEM_antiPrefStops);
%             SEMhigh = [meanantiPrefStops(keepIndex) + SEM_antiPrefStops(keepIndex)]; 
%             SEMlow = [meanantiPrefStops(keepIndex) - SEM_antiPrefStops(keepIndex)];
%             ax2 = subplot(2,1,2);
%             plot(windowTime, meanantiPrefStops,'k','LineWidth',1.5)
%             patch([windowTime(keepIndex) fliplr(windowTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%             title(ax2,'antiPrefStops')
%             xline(ax2,0)
%             box off
%             %ylim(ax2,[-63,-56])
%             linkaxes([ax1,ax2],'x')
%             saveas(gcf,fullfile(fly,'StopTransition_sum.fig'))
% 
%             figure(); 
%             set(gcf,'color','w')
%             set(gcf,'renderer','painters')
%             sgtitle(temp(breaks(end-2) + 1:end))
%             keepIndex = ~isnan(SEM_prefStarts);
%             SEMhigh = [meanPrefStarts(keepIndex) + SEM_prefStarts(keepIndex)]; 
%             SEMlow = [meanPrefStarts(keepIndex) - SEM_prefStarts(keepIndex)];
%             bx1 = subplot(2,1,1);
%             plot(windowTime,meanPrefStarts,'k','LineWidth',1.5)
%             patch([windowTime(keepIndex) fliplr(windowTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%             title(bx1,'prefStarts')
%             xline(bx1,0)
%             box off
%             %ylim(bx1,[-63,-56])
%             keepIndex = ~isnan(SEM_antiPrefStarts);
%             SEMhigh = [meanantiPrefStarts(keepIndex) + SEM_antiPrefStarts(keepIndex)]; 
%             SEMlow = [meanantiPrefStarts(keepIndex) - SEM_antiPrefStarts(keepIndex)];
%             bx2 = subplot(2,1,2);
%             plot(windowTime,meanantiPrefStarts,'k','LineWidth',1.5)
%             patch([windowTime(keepIndex) fliplr(windowTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%             title(bx2,'antiPrefStarts')
%             xline(bx2,0)
%             box off
%             linkaxes([bx1,bx2],'x')
%             
%             saveas(gcf,fullfile(fly,'StartTransition_sum.fig'))
            
            flyPrefStops(t,:) = mean(prefStops{t},1);
            flyantiPrefStops(t,:) = mean(antiPrefStops{t},1);
            flyPrefStarts(t,:) = mean(prefStarts{t},1);
            flyantiPrefStarts(t,:) = mean(antiPrefStarts{t},1);

end

            aveVmprefStops = mean(flyPrefStops(:,1));
            aveVmantiPrefStops = mean(flyantiPrefStops(:,1));
            aveVmprefStarts = mean(flyPrefStarts(:,1));
            aveVmantiPrefStarts = mean(flyantiPrefStarts(:,1));
            
%             baselineDiff = flyPrefStops(:,1) - aveVmprefStops;
%             flyPrefStops = flyPrefStops - baselineDiff;
%             
%             baselineDiff = flyantiPrefStops(:,1) - aveVmantiPrefStops;
%             flyantiPrefStops = flyantiPrefStops - baselineDiff;
%             
%             baselineDiff = flyPrefStarts(:,1) - aveVmprefStarts;
%             flyPrefStarts = flyPrefStarts - baselineDiff;
%             
%             baselineDiff = flyantiPrefStarts(:,1) - aveVmantiPrefStarts;
%             flyantiPrefStarts = flyantiPrefStarts - baselineDiff;
%             
% 
%             SEM_prefStops = std(flyPrefStops,[],1,'omitnan') / sqrt(size(flyPrefStops,1));
%             SEM_antiPrefStops = std(flyantiPrefStops,[],1,'omitnan') / sqrt(size(flyantiPrefStops,1));
%             SEM_prefStarts = std(flyPrefStarts,[],1,'omitnan') / sqrt(size(flyPrefStarts,1));
%             SEM_antiPrefStarts = std(flyantiPrefStarts,[],1,'omitnan') / sqrt(size(flyantiPrefStarts,1));
%     
%             meanPrefStops = mean(flyPrefStops,1);
%             meanantiPrefStops = mean(flyantiPrefStops,1);
%             meanPrefStarts = mean(flyPrefStarts,1);
%             meanantiPrefStarts = mean(flyantiPrefStarts,1);
%             
%             figure();
%             set(gcf,'color','w')
%             set(gcf,'renderer','painters')
%             sgtitle('all flies')
%             keepIndex = ~isnan(SEM_prefStops);
%             SEMhigh = [meanPrefStops(keepIndex) + SEM_prefStops(keepIndex)]; 
%             SEMlow = [meanPrefStops(keepIndex) - SEM_prefStops(keepIndex)];
%             ax1 = subplot(2,1,1);
%             plot(windowTime, meanPrefStops,'k','LineWidth',1.5)
%             patch([windowTime(keepIndex) fliplr(windowTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%             title(ax1,'prefStops')
%             xline(ax1,0)
%             box off
%             keepIndex = ~isnan(SEM_antiPrefStops);
%             SEMhigh = [meanantiPrefStops(keepIndex) + SEM_antiPrefStops(keepIndex)]; 
%             SEMlow = [meanantiPrefStops(keepIndex) - SEM_antiPrefStops(keepIndex)];
%             ax2 = subplot(2,1,2);
%             plot(windowTime, meanantiPrefStops,'k','LineWidth',1.5)
%             patch([windowTime(keepIndex) fliplr(windowTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%             title(ax2,'antiPrefStops')
%             xline(ax2,0)
%             box off
%             %ylim(ax2,[-63,-56])
%             linkaxes([ax1,ax2],'x')
%             
%             saveas(gcf,fullfile('Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL3\goodFlies','StopTransition_allFlies.fig'))
% 
%             figure(); 
%             set(gcf,'color','w')
%             set(gcf,'renderer','painters')
%             sgtitle('all flies')
%             keepIndex = ~isnan(SEM_prefStarts);
%             SEMhigh = [meanPrefStarts(keepIndex) + SEM_prefStarts(keepIndex)]; 
%             SEMlow = [meanPrefStarts(keepIndex) - SEM_prefStarts(keepIndex)];
%             bx1 = subplot(2,1,1);
%             plot(windowTime,meanPrefStarts,'k','LineWidth',1.5)
%             patch([windowTime(keepIndex) fliplr(windowTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%             title(bx1,'prefStarts')
%             xline(bx1,0)
%             box off
%             %ylim(bx1,[-63,-56])
%             keepIndex = ~isnan(SEM_antiPrefStarts);
%             SEMhigh = [meanantiPrefStarts(keepIndex) + SEM_antiPrefStarts(keepIndex)]; 
%             SEMlow = [meanantiPrefStarts(keepIndex) - SEM_antiPrefStarts(keepIndex)];
%             bx2 = subplot(2,1,2);
%             plot(windowTime,meanantiPrefStarts,'k','LineWidth',1.5)
%             patch([windowTime(keepIndex) fliplr(windowTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%             title(bx2,'antiPrefStarts')
%             xline(bx2,0)
%             box off
%             linkaxes([bx1,bx2],'x')
%             
%             saveas(gcf,fullfile('Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL3\goodFlies','StartTransition_allFlies.fig'))
%%              cells = unique(jumpSum.expID);
    
    nCol = length(antiPrefStarts); 
    colormap = cbrewer2('RdBu',nCol+1);  
    colormap(1,:) = [];
    colormap_new = zeros(size(colormap)); 
    colormap_new(1,:) = colormap(1,:); 
    colormap_new(2,:) = colormap(1,:);
    colormap_new(3,:) = colormap(4,:); 
    colormap_new(4,:) = colormap(4,:); 
    colormap_new(5,:) = colormap(4,:); 
    colormap_new(6,:) = colormap(7,:); 
    colormap_new(7,:) = colormap(7,:); 
    figure();    
    hold on
    for fly = 1:nCol
        hold on
        set(gcf,'color','w')
        set(gcf,'renderer','painters')
        %sgtitle('all flies')
        %ax1 = subplot(2,1,1);
        plot(windowTime, flyPrefStops(fly,:),'color',colormap_new(fly,:),'LineWidth',1.5)
        title('prefStops')
        xline(0)
        box off
    end
    figure();    
    hold on
    for fly = 1:nCol
        %ax2 = subplot(2,1,2);
        set(gcf,'color','w')
        set(gcf,'renderer','painters')
        plot(windowTime, flyantiPrefStops(fly,:),'color',colormap_new(fly,:),'LineWidth',1.5)
        title('antiPrefStops')
        xline(0)
        box off
        %ylim(ax2,[-63,-56])
        %linkaxes([ax1,ax2],'x')
    end
            
    figure();    
    hold on
    for fly = 1:nCol 
        set(gcf,'color','w')
        set(gcf,'renderer','painters')
        %sgtitle('all flies')
        %bx1 = subplot(2,1,1);
        plot(windowTime,flyPrefStarts(fly,:),'color',colormap_new(fly,:),'LineWidth',1.5)
        title('prefStarts')
        xline(0)
        box off
    end
    figure();    
    hold on
    for fly = 1:nCol
        set(gcf,'color','w')
        set(gcf,'renderer','painters')
        %ylim(bx1,[-63,-56])
        %bx2 = subplot(2,1,2);
        plot(windowTime,flyantiPrefStarts(fly,:),'color',colormap_new(fly,:),'LineWidth',1.5)
        title('antiPrefStarts')
        xline(0)
        box off
        %linkaxes([bx1,bx2],'x')
    end


end