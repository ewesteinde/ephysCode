function activityVSbehaviour_ionto_mean(folders,behaviour_data,trialData, data, savePlot)       
    edges_vy = [0:5:150];
    edges_vf = [-4:1:10];
    edges_angle = [-180:10:180]; 
    
    threshold=3;
    trial_count = 1;
    
    trial_vy = [];
    trial_angle = [];
    trial_vf = []; 
    speed_sum = [];
    
for f = 1:size(folders,1)
    folder = folders{f};
        try
%% Remove idx where the fly isn't moving
        total_mov_mm = abs(behaviour_data.vel_for + abs(behaviour_data.vel_side) + abs(deg2rad(behaviour_data.vel_yaw))*4.5);
        no0vel_idx = find(total_mov_mm > threshold);
        vf = behaviour_data.vel_for;
        vf = vf(no0vel_idx); 
        vs = behaviour_data.vel_side;
        vs = abs(vs(no0vel_idx)); 
        vy = behaviour_data.vel_yaw;
        vy = abs(vy(no0vel_idx)); 
        angle = behaviour_data.angle; 
        angle = angle(no0vel_idx); 
        %angle = wrapTo180(wrapTo180(angle)-wrapTo180(goal); 

    %%

        sum_mean{1} = zeros(length(edges_vy)-1,1);
        sum_mean{2} = zeros(length(edges_angle)-1,1);
        sum_mean{3} = zeros(length(edges_vf)-1,1);

        if strcmp(data,'fRate')
            activity = trialData.fRate_sec; 
        else
            activity = trialData.smoothVm;
        end
                activity = activity(no0vel_idx);

                % vf 
                behaviour = vf; 
                [zscore, centers_vf, ~] = binData(activity, behaviour, edges_vf);
                sum_mean{3} = zscore;
                
                % vy 
                 behaviour = vy; 
                [zscore, centers_vy, ~] = binData(activity, behaviour, edges_vy);
                sum_mean{1} = zscore; 
                

                % angle
                [zscore, centers_angle, ~] = binData(activity, angle, edges_angle);
                sum_mean{2}= zscore; 
        
        
        [~, vy_vf_activity, heatvy_centers, heatvf_centers] = create_activity_heatmap(vy, vf, activity, edges_vy, edges_vf);
        [~, angle_vf_activity, heatangle_centers, heatvf_centers] = create_activity_heatmap(angle, vf, activity, edges_angle, edges_vf);
        trial_vf_vy(trial_count,:,:) = vy_vf_activity;
        trial_vf_angle(trial_count,:,:) = angle_vf_activity;
        
        trial_vf(trial_count,:) = sum_mean{3};
        trial_vy(trial_count,:) = sum_mean{1};
        trial_angle(trial_count,:) = sum_mean{2};

        trial_count = trial_count + 1; 
        
        catch
            disp(['folder ',folder,' failed'])
        end
end

%% heatmaps    
            if ~strcmp(data,'fRate')
                trial_vf_vy(trial_vf_vy == 0) = nan; 
                trial_vf_angle(trial_vf_angle == 0) = nan; 

                trial_vf(trial_vf == 0) = nan; 
                trial_vy(trial_vy == 0) = nan; 
                trial_angle(trial_angle == 0) = nan; 
            end


            aveTrial_vy_vf = squeeze(mean(trial_vf_vy,1,'omitnan'));
            aveTrial_angle_vf = squeeze(mean(trial_vf_angle,1,'omitnan'));
            
            aveTrial_vy_vf(aveTrial_vy_vf == 0) = nan;
            ncol = length(unique(aveTrial_vy_vf));
            color = flipud(cbrewer2('RdYlBu', ncol));
            figure();
            set(gcf,'color','w')
            set(gcf,'renderer','painters')
            subplot(2,1,1);
            s = pcolor(aveTrial_vy_vf); 
            colormap(color)
            ylabel('Vf mm/s')
            xt = linspace(1,numel(heatvy_centers),7);                            
            xtlbl = linspace(heatvy_centers(xt(1)), heatvy_centers(xt(end)), 7);   
            set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
            xlabel('Vy deg/s')
            colorbar
            yt = linspace(1,numel(heatvf_centers),5); 
            ytlbl = linspace(heatvf_centers(yt(1)), heatvf_centers(yt(end)), 5);
            set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
            set(s, 'EdgeColor', 'none');
            set(gca,'color','none')
            box off
            
            aveTrial_angle_vf(aveTrial_angle_vf == 0) = nan;
            ncol = length(unique(aveTrial_angle_vf));
            color = flipud(cbrewer2('RdYlBu', ncol));
            subplot(2,1,2);
            s = pcolor(aveTrial_angle_vf); 
            colormap(color)
            ylabel('Vf mm/s')
            xt = linspace(1,numel(heatangle_centers),5);                            
            xtlbl = linspace(heatangle_centers(xt(1)), heatangle_centers(xt(end)), 5);   
            set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
            xlabel('cue pos')
            colorbar
            yt = linspace(1,numel(heatvf_centers),5); 
            ytlbl = linspace(heatvf_centers(yt(1)), heatvf_centers(yt(end)), 5);
            set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
            set(s, 'EdgeColor', 'none');
            set(gca,'color','none')
            box off
            
        if savePlot
            saveas(gcf, fullfile(folder,'figures',[data,'_activityvsBehaviour_heatmap.fig']));
            saveas(gcf, fullfile(folder,'figures',[data,'_activityvsBehaviour_heatmap.svg']));
        end
             
%% lineplots    
%             allTrial_vf = var(trial_vf,[],3,'omitnan');
%             allTrial_vy = var(trial_vy,[],3,'omitnan');
%             allTrial_angle = var(trial_angle,[],3,'omitnan');
            
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
            ylabel(data)
            
            
            keepIndex = ~isnan(SEM_angle);
            SEMhigh = [aveTrial_angle(keepIndex) + SEM_angle(keepIndex)]; 
            SEMlow = [aveTrial_angle(keepIndex) - SEM_angle(keepIndex)];
            ax2 = subplot(3,1,2);
            hold on
            plot(centers_angle, aveTrial_angle,'b')
            patch([centers_angle(keepIndex) fliplr(centers_angle(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
            ax1.YAxis.Exponent = 0;
            ax2.YAxis.Exponent = 0;
            xlabel('cue pos')
            ylabel(data)

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
            ylabel(data)
            ax3.YAxis.Exponent = 0;

        if savePlot
            saveas(gcf, fullfile(folder,'figures',[data,'_activityvsBehaviour_lineplots.fig']));
            saveas(gcf, fullfile(folder,'figures',[data,'_activityvsBehaviour_lineplots.svg']));
        end
      
end