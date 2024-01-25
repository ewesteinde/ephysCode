function activityVSbehaviour_ephysFlyPrefHeadSummary(summaryArray,meno,fRate)       
    edges_vy = [-200:10:200];
    edges_vf = [-4:1:10];
    edges_vs = [-8:1:8];
    edges_angle = [-180:30:180]; 
    
    threshold=2;
    trial_count = 1;
    
    trial_vy = [];
    trial_angle = [];
    trial_vf = []; 
    speed_sum = [];
    
    activity_values = [];

    summaryArray = summaryArray(~ismembertol(summaryArray.rho,1,10^-10),:); % gets rid of trials w/ no heading change --> indicates problem
    summaryArray = summaryArray(summaryArray.timeMov > 10,:); % only look at trials where fly's vel was above threshold for at least 5 seconds

%     if meno
%         summaryArray = summaryArray(summaryArray.rho > 0.5,:);
%     else
%         summaryArray = summaryArray(summaryArray.rho < 0.5,:);
%     end

for f = 1:size(summaryArray,1)
    folder = summaryArray.Folder(f); 
    trialIdx = regexp(folder,'_trial');
    if ~isempty(trialIdx)
        [~,cellIdx] = regexp(folder,'cell_');
        folder_new = folder{1}(1:cellIdx + 1);
        flies(f) = string(folder_new);
    else
        flies(f) = string(folder); 
    end    
end

flies = flies';

    for chunk = 1:size(summaryArray,1)
        try
        folder = table2array(summaryArray(chunk,1)); 
        if strcmp(folder(end),'.')
            folder = folder(1:end-2); 
        end

        processedDir = fullfile(folder,'processedData');
        load(fullfile(processedDir,'pro_behaviourData.mat'),'pro_behaviourData')
        load(fullfile(processedDir,'pro_trialData.mat'),'pro_trialData')
    
        bData = pro_behaviourData{summaryArray.numTrial(chunk)};
        tData = pro_trialData{summaryArray.numTrial(chunk)}; 

% Remove idx where the fly isn't moving
        total_mov_mm = abs(bData.vel_for(summaryArray.Indices{chunk}) + abs(bData.vel_side(summaryArray.Indices{chunk})) + abs(deg2rad(bData.vel_yaw(summaryArray.Indices{chunk}))*4.5));
        no0vel_idx = find(total_mov_mm > threshold);
        vf = bData.vel_for(summaryArray.Indices{chunk});
        vf = vf(no0vel_idx); 
        vs = bData.vel_side(summaryArray.Indices{chunk});
        vs = vs(no0vel_idx); 
        vy = bData.vel_yaw(summaryArray.Indices{chunk});
        vy = vy(no0vel_idx); 
        speed = sqrt(vf.^2 + vs.^2); 
        angle = bData.angle(summaryArray.Indices{chunk}); 
        angle = angle(no0vel_idx); 
        %angle = wrapTo180(wrapTo180(angle)-wrapTo180(rad2deg(summaryArray.Goal(trial)))); 
        speed_sum = [speed_sum,speed'];


        sum_mean{2} = zeros(length(edges_angle)-1,1);
        sum_mean{3} = zeros(length(edges_vf)-1,1);
        sum_mean{1} = zeros(length(edges_vy)-1,1);

        if fRate    
            activity = tData.fRate_sec;
            activity = activity(summaryArray.Indices{chunk}); 
            activity = activity(no0vel_idx);
        else
            try
                activity = tData.smooth_Vm;
            catch
                activity = tData.smoothVm;
            end
            activity = activity(summaryArray.Indices{chunk}); 
            activity = activity(no0vel_idx);
        end
        
        activity = (activity - min(activity))/(max(activity)- min(activity));

        %vf 
        behaviour = vf; 
        [zscore, centers_vf, ~] = binData(activity, behaviour, edges_vf);
        sum_mean{3} = zscore;

       % vy 
         behaviour = vy; 
        [zscore, centers_vy, ~] = binData(activity, behaviour, edges_vy);
        sum_mean{1} = zscore; 

        %angle
        [zscore, centers_angle, ~] = binData(activity, angle, edges_angle);
        sum_mean{2} = zscore; 
        
%         [N_vy, vy_vf_activity, heatvy_centers, heatvf_centers] = create_activity_heatmap(vy, vf, activity, edges_vy, edges_vf);
%         [N_Cue, angle_vf_activity, heatangle_centers, heatvf_centers] = create_activity_heatmap(angle, vf, activity, edges_angle, edges_vf);
%         trial_vf_vy(trial_count,:,:) = vy_vf_activity;
%         trial_vf_angle(trial_count,:,:) = angle_vf_activity;
%         N_totalCue = N_totalCue + N_Cue;
%         N_totalvy = N_totalvy + N_vy;

        prefHead_smooth = smoothdata(sum_mean{2},'gaussian',5);
        
        prefHead = centers_angle(prefHead_smooth == max(prefHead_smooth));
        
%         bump_params = fit_sinusoid(prefHead_smooth, [0,2*pi], 1);
%         prefHead = wrapTo180((rad2deg(wrapToPi(bump_params.pos_rad + deg2rad(15)))) - 180);
% %         figure(1);clf;
%         plot(centers_angle,sum_mean{2})
%         hold on
%         plot(centers_angle, prefHead_smooth)
%         xline(prefHead)
        
        if fRate
            if sum(~isnan(sum_mean{2})) < length(sum_mean{2}) - 2
                FRprefHead(chunk) = nan; 
                FRprefHead_amp(chunk) = nan;
            else
                FRprefHead(chunk) = prefHead; 
                FRprefHead_amp(chunk) = max(prefHead_smooth) - min(prefHead_smooth); %bump_params.amp;
            end
        else
            if sum(~isnan(sum_mean{2})) < length(sum_mean{2})  - 2
                VMprefHead(chunk) = nan; 
                VMprefHead_amp(chunk) = nan;
            else
                VMprefHead(chunk) = prefHead; 
                VMprefHead_amp(chunk) = max(prefHead_smooth) - min(prefHead_smooth); %bump_params.amp;
            end
        end
        
        trial_vf(trial_count,:) = sum_mean{3};
        trial_vy(trial_count,:) = sum_mean{1};
        trial_angle(trial_count,:) = sum_mean{2};

        trial_count = trial_count + 1; 
        
       % pause(0.2)
        catch
            disp(['folder ',folder,' failed'])
%             if fRate
%                 FRprefHead(chunk) = nan; 
%                 FRprefHead_amp(chunk) = nan;
%             else
%                 VMprefHead(chunk) = nan; 
%                 VMprefHead_amp(chunk) = nan;
%             end
%             trial_vf(chunk,:) = nan;
%             trial_vy(chunk,:) = nan;
%             trial_angle(chunk,:) = nan;
        end
    end  


% heatmaps    
            
%             trial_vf_vy(trial_vf_vy == 0) = nan; 
%             trial_vf_angle(trial_vf_angle == 0) = nan; 
% 
%             trial_vf(trial_vf == 0) = nan; 
%             trial_vy(trial_vy == 0) = nan; 
%             trial_angle(trial_angle == 0) = nan; 
% 
% 
%             aveTrial_vy_vf = squeeze(mean(trial_vf_vy,1,'omitnan'));
%             aveTrial_angle_vf = squeeze(mean(trial_vf_angle,1,'omitnan'));
%             greyBins = N_totalvy < 100 ;
%             aveTrial_vy_vf(aveTrial_vy_vf == 0) = nan;
%             aveTrial_vy_vf(greyBins) = nan; 
%             ncol = length(unique(aveTrial_vy_vf));
%             color = flipud(cbrewer2('RdYlBu', ncol));
%             figure();
%             set(gcf,'color','w')
%             set(gcf,'renderer','painters')
%             subplot(2,1,1);
%             s = pcolor(aveTrial_vy_vf);
%             colormap(color)
%             hold on
% 
%             greyBins(greyBins == 0) = nan; 
%             pcolor(greyBins)
%             cmap = ['none',0.5,0.5,0.5;]
%             
%             
%             ylabel('Vf mm/s')
%             xt = linspace(1,numel(heatvy_centers),7);                            
%             xtlbl = linspace(heatvy_centers(xt(1)), heatvy_centers(xt(end)), 7);   
%             set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
%             xlabel('Vy deg/s')
%             colorbar
%             yt = linspace(1,numel(heatvf_centers),5); 
%             ytlbl = linspace(heatvf_centers(yt(1)), heatvf_centers(yt(end)), 5);
%             set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
%             set(s, 'EdgeColor', 'none');
%             set(gca,'color','none')
%             box off
%             
%             greyBins = N_totalCue < 100 ;
%             aveTrial_angle_vf(aveTrial_angle_vf == 0) = nan;
%             aveTrial_angle_vf(greyBins) = nan;
%             ncol = length(unique(aveTrial_angle_vf));
%             color = flipud(cbrewer2('RdYlBu', ncol));
%             subplot(2,1,2);
%             s = pcolor(aveTrial_angle_vf); 
%             colormap(color)
%             ylabel('Vf mm/s')
%             xt = linspace(1,numel(heatangle_centers),5);                            
%             xtlbl = linspace(heatangle_centers(xt(1)), heatangle_centers(xt(end)), 5);   
%             set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
%             xlabel('cue pos rel to goal')
%             colorbar
%             yt = linspace(1,numel(heatvf_centers),5); 
%             ytlbl = linspace(heatvf_centers(yt(1)), heatvf_centers(yt(end)), 5);
%             set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
%             set(s, 'EdgeColor', 'none');
%             set(gca,'color','none')
%             box off
            
             
% lineplots    
            
            if fRate
                summaryArray.FRprefHead = FRprefHead'; 
                summaryArray.FRprefHead_amp =  FRprefHead_amp';
            else
                summaryArray.VMprefHead = VMprefHead'; 
                summaryArray.VMprefHead_amp =  VMprefHead_amp';
            end
            
            allFlies = unique(flies);
            if fRate
                goalHeadDiff = rad2deg(angdiff(summaryArray.Goal,deg2rad(summaryArray.FRprefHead)));
                summaryArray.FRgoalHeadDiff = goalHeadDiff; 
                figure(); 
                pointsize = 20; 
                scatter(summaryArray.FRgoalHeadDiff,summaryArray.FRprefHead_amp,pointsize,summaryArray.rho,'filled')
                
                for i = 1:length(allFlies)
                    fly = allFlies(i); 
                    flyData = summaryArray(flies == fly,:); 


                    figure(); 
                    pointsize = 10; 
                    scatter(flyData.FRgoalHeadDiff,flyData.FRprefHead_amp,pointsize,flyData.rho,'filled')
                    xlabel('|goal - prefHead|')
                    ylabel('prefHead amp')
                    colorbar
                    figure();
                    scatter(flyData.rho, flyData.FRprefHead_amp)
                    xlabel('rho')
                    ylabel('prefHead amp')
                end
            else
                goalHeadDiff = rad2deg(angdiff(summaryArray.Goal,deg2rad(summaryArray.VMprefHead)));
                summaryArray.VMgoalHeadDiff = goalHeadDiff; 
                figure(); 
                pointsize = 30; 
                scatter(summaryArray.VMgoalHeadDiff,summaryArray.VMprefHead_amp,pointsize,summaryArray.rho,'filled')
                

                for i = 1:length(allFlies)
                    fly = allFlies(i); 
                    flyData = summaryArray(flies == fly,:); 


                    figure(); 
                    pointsize = 10; 
                    scatter(flyData.VMgoalHeadDiff,flyData.VMprefHead_amp,pointsize,flyData.rho,'filled')
                    xlabel('|goal - prefHead|')
                    ylabel('prefHead amp')
                    colorbar
                    figure();
                    scatter(flyData.rho, flyData.VMprefHead_amp)
                    xlabel('rho')
                    ylabel('prefHead amp')
                end
            end

            allTrial_vf = var(trial_vf,[],3,'omitnan');
            allTrial_vy = var(trial_vy,[],3,'omitnan');
            allTrial_angle = var(trial_angle,[],3,'omitnan');
            
            SEM_vf = std(trial_vf,[],1,'omitnan') / sqrt(size(trial_vf,1));
            SEM_vy = std(trial_vy,[],1,'omitnan') / sqrt(size(trial_vy,1));
            SEM_angle = std(trial_angle,[],1,'omitnan') / sqrt(size(trial_angle,1));
            
            aveTrial_vf = mean(trial_vf,1,'omitnan');
            aveTrial_vy = mean(trial_vy,1,'omitnan');
            aveTrial_angle = mean(trial_angle,1,'omitnan');
            
            
            figure();
            set(gcf,'color','w')
            set(gcf,'Renderer','painters')
            keepIndex = ~isnan(SEM_vy);
            SEMhigh = [aveTrial_vy(keepIndex) + SEM_vy(keepIndex)]; 
            SEMlow = [aveTrial_vy(keepIndex) - SEM_vy(keepIndex)];
            ax1 = subplot(3,1,1);
            hold on
            plot(centers_vy, aveTrial_vy,'b')
            patch([centers_vy(keepIndex) fliplr(centers_vy(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
            xlabel('vy (deg/s)')
            ylabel('bump amplitude')
            
            keepIndex = ~isnan(SEM_angle);
            SEMhigh = [aveTrial_angle(keepIndex) + SEM_angle(keepIndex)]; 
            SEMlow = [aveTrial_angle(keepIndex) - SEM_angle(keepIndex)];
            ax2 = subplot(3,1,2);
            hold on
            plot(centers_angle, aveTrial_angle,'b')
            patch([centers_angle(keepIndex) fliplr(centers_angle(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
            ax1.YAxis.Exponent = 0;
            ax2.YAxis.Exponent = 0;
            xlabel('cue pos rel goal')
            ylabel('bump amplitude')

            set(gcf,'color','w')
            set(gcf,'Renderer','painters')
            keepIndex = ~isnan(SEM_vf);
            SEMhigh = [aveTrial_vf(keepIndex) + SEM_vf(keepIndex)]; 
            SEMlow = [aveTrial_vf(keepIndex) - SEM_vf(keepIndex)];
            ax3 = subplot(3,1,3);
            hold on
            plot(centers_vf, aveTrial_vf,'b')
            patch([centers_vf(keepIndex) fliplr(centers_vf(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
            xlabel('vf (mm/s)')
            ylabel('bump amplitude')
            ax3.YAxis.Exponent = 0;

    if savePlots == 1
        saveas(z, fullfile(lineplotDir,[expID,'_',num2str(nTrial),'_zScore_behaviour_no0vel.fig']));
    end
end