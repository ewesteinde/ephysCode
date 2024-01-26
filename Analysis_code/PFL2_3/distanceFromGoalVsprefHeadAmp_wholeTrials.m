folders = get_folders_ephys(rootDir);

savePlots = 1; 
count = 1; 
rho = []; 
vecStr = []; 
prefHead_all = []; 
prefHead_amp = []; 
prefHead_smooth = [];
goal = []; 
goalHeadDiff = []; 
for f = 1:size(folders,1)
    try
        folder = folders(f).folder; 
        if strcmp(folder(end),'.')
            folder = folder(1:end-2); 
        end

        processedDir = fullfile(folder,'processedData');
        load(fullfile(processedDir,'pro_behaviourData.mat'))
        load(fullfile(processedDir,'pro_trialData.mat'))

        numTrials = length(pro_behaviourData);

        %% visualize individual trials

        for t = 1:numTrials
            bData = pro_behaviourData{t};
            tData = pro_trialData{t};
            [rho, theta] = CalculateAverageHeading_ephys(bData,1.5, 'all');
            goal(count) = -theta; % - to convert to 'heading' to match model
            vecStr(count) = rho; 
            
             total_mov_mm = abs(bData.vel_for + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw)*4.5));
             no0vel_idx = find(total_mov_mm > 2);
             angle = -bData.angle(no0vel_idx); % - to convert to 'heading' to match model
            
            if fRate    
                activity = tData.fRate_sec;
                activity = activity(no0vel_idx);
                activity = (activity - min(activity))/(max(activity)- min(activity));
            else
                try
                    activity = tData.smooth_Vm;
                catch
                    activity = tData.smoothVm;
                end
                activity = activity(no0vel_idx);
                activity = (activity - min(activity))/(max(activity)- min(activity));
            end
            edges_angle = [-180:30:180]; 
            [activity_angle, centers_angle, ~] = binData(activity, angle, edges_angle);
            
            prefHead_smooth = smoothdata(activity_angle,'gaussian',5);
            %prefHead = centers_angle(prefHead_smooth == max(prefHead_smooth));
            
            bump_params = fit_sinusoid(prefHead_smooth, [0,2*pi], 1);
            prefHead = wrapTo180((rad2deg(wrapToPi(bump_params.pos_rad + deg2rad(15)))) - 180);

%             figure(1);clf;
%             plot(centers_angle,activity_angle)
%             hold on
%             plot(centers_angle, prefHead_smooth)
%             xline(prefHead)
            
            if sum(~isnan(activity_angle)) < length(activity_angle) - 2 || bump_params.adj_rs < 0.5
                prefHead_all(count) = nan; 
                prefHead_amp(count) = nan;
            else
                prefHead_all(count) = prefHead; 
                prefHead_amp(count) = bump_params.amp;%max(prefHead_smooth) - min(prefHead_smooth);
            end

            count = count + 1; 
        end
    catch
        disp([folder,' failed'])
                prefHead_all(count) = nan; 
                prefHead_amp(count) = nan;
        count = count + 1; 
    end  
end
goalHeadDiff = rad2deg(angdiff(goal,deg2rad(prefHead_all)));
% [N_Vm, heatmap_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(abs(goalHeadDiff)', vecStr', prefHead_amp', [0:30:180], [0:0.1:1]);
% 
%             ncol = length(unique(heatmap_Vm));
%             color = flipud(cbrewer2('RdYlBu', ncol));
%             figure();
%             set(gcf,'color','w')
%             set(gcf,'renderer','painters')
%             s = pcolor(heatmap_Vm); 
%             colormap(color)
%             ylabel('vector str')
%             xt = linspace(1,numel(x_centers_Vm),7);                            
%             xtlbl = linspace(x_centers_Vm(xt(1)), x_centers_Vm(xt(end)), 7);   
%             set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
%             xlabel('goal - prefHead')
%             colorbar
%             yt = linspace(1,numel(y_centers_Vm),5); 
%             ytlbl = linspace(y_centers_Vm(yt(1)), y_centers_Vm(yt(end)), 5);
%             set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
%             set(s, 'EdgeColor', 'none');
%             set(gca,'color','none')
%             box off
            
            

figure(); 
set(gcf,'color','w')
set(gcf,'renderer','painters')
scatter(goalHeadDiff',prefHead_amp',20,vecStr','filled')
a = colorbar;
ylabel(a,'vector str')
xlabel('goal - prefHead')
ylabel('prefHead amp')
