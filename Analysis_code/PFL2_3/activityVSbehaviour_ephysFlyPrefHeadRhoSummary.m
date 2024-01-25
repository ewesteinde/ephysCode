function activityVSbehaviour_ephysFlyPrefHeadRhoSummary(rootDir,vm)       
    edges_vy = [-200:10:200];
    edges_vf = [-4:1:10];
    edges_vs = [-8:1:8];
    edges_angle = [-180:30:180]; 
    
    threshold=1.5;
    trial_count = 1;
    
    trial_vy = [];
    trial_angle = [];
    trial_vf = []; 
    speed_sum = [];
    
    activity_values = [];
    
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
count = 1; 

for i = 1:length(allFlies)
    
    fly = allFlies(i);
    
    fly_folders = folders(flies == fly);
try
    for f = 1:size(fly_folders,1)
        Meno_angle = [];
        nMeno_angle = [];
        Meno_activity = [];
        nMeno_activity = [];
        prefHead = [];   
        meno_idx = [];
        nMeno_idx = [];
    
        folder = fly_folders(f).folder; 
        if strcmp(folder(end),'.')
            folder = folder(1:end-2); 
        end
        
        BreakIdx = regexp(folder, '\');
        cellSide = folder(BreakIdx(7)-1);
        clear pro_behaviourData
        clear processed_behaviourData
        
        clear pro_trialData
        clear processed_trialData

        processedDir = fullfile(folder,'processedData');
        try
            load(fullfile(processedDir,'pro_behaviourData.mat'))
            load(fullfile(processedDir,'pro_trialData.mat'))
        catch

            load(fullfile(folder,'pro_behaviourData.mat'))
            load(fullfile(folder,'pro_trialData.mat'))
        end
        load(fullfile(folder,'MenoDataFly_trial_1.mat'))
        load(fullfile(folder,'nMenoDataFly_trial_1.mat'))
    
        try
            bData = pro_behaviourData{1};
            tData = pro_trialData{1};
        catch
            bData = processed_behaviourData{1};
            tData = processed_trialData{1};
        end

        prefHead(f) = calcPrefHead(bData, tData);
        if vm           
            try
                activity = tData.smooth_Vm;
            catch
                activity = tData.smoothVm;
            end
        else        
            activity = tData.fRate_sec;
        end

        for c = 1:size(MenoDataFly,1)
        % Remove idx where the fly isn't moving
            Meno_total_mov_mm = abs(bData.vel_for(MenoDataFly.Indices{c}) + abs(bData.vel_side(MenoDataFly.Indices{c})) + abs(deg2rad(bData.vel_yaw(MenoDataFly.Indices{c}))*4.5));
            Meno_no0vel_idx = find(Meno_total_mov_mm > threshold);
            chunk_angle = bData.angle(MenoDataFly.Indices{c}); 
            meno_idx = [meno_idx; MenoDataFly.Indices{c}'];
            Meno_angle = [Meno_angle; chunk_angle(Meno_no0vel_idx)]; 
            
            chunk_activity = activity(MenoDataFly.Indices{c}); 
            Meno_activity = [Meno_activity; chunk_activity(Meno_no0vel_idx)];
        end
        
        for c = 1:size(nMenoDataFly,1)
            nMeno_total_mov_mm = abs(bData.vel_for(nMenoDataFly.Indices{c}) + abs(bData.vel_side(nMenoDataFly.Indices{c})) + abs(deg2rad(bData.vel_yaw(nMenoDataFly.Indices{c}))*4.5));
            nMeno_no0vel_idx = find(nMeno_total_mov_mm > threshold);
            chunk_angle = bData.angle(nMenoDataFly.Indices{c}); 
            nMeno_angle = [nMeno_angle; chunk_angle(nMeno_no0vel_idx)]; 
            nMeno_idx = [nMeno_idx; nMenoDataFly.Indices{c}'];
            chunk_activity = activity(nMenoDataFly.Indices{c});
            nMeno_activity = [nMeno_activity; chunk_activity(nMeno_no0vel_idx)];
        end  
        %activity = (activity - min(activity))/(max(activity)- min(activity));
        
        [Meno_prefHead, centers_angle, ~] = binData(Meno_activity, Meno_angle, edges_angle);
        
        [nMeno_prefHead, ~, ~] = binData(nMeno_activity, nMeno_angle, edges_angle);
        
        allGoals = repelem(MenoDataFly.Goal,round(MenoDataFly.timeMov)); 
        if size(allGoals,2) > size(allGoals,1)
            allGoals = allGoals';
        end

        ave_goal = -(circ_mean(allGoals)); % switch from cue pos to HD
        avePrefHead = -mean(prefHead,'omitnan');

        goalDiff = rad2deg(angdiff(ave_goal,deg2rad(avePrefHead))); 

        figure(11);clf;
        set(gcf,'color','w','renderer','painters')
        a1 = subplot(1,1,1);
        hold on
        plot(-centers_angle,Meno_prefHead,'k')
        plot(-centers_angle,nMeno_prefHead,'color',[0.75,0.75,0.75])
        xlabel('HD')
        ylabel('Vm')

        if regexp(rootDir,'PFL3')
            title([cellSide,' goalDiff: ',num2str(goalDiff)])
        else
            title('goalDiff: ',num2str(goalDiff))
        end
        box off
        xlim(a1,[-180,180])
        maxMeno = max(Meno_prefHead);
        minMeno = min(Meno_prefHead);
    end
    
catch
    disp([folder,' failed'])
end
        %Meno_prefHead_smooth = smoothdata(Meno_binAngle,'gaussian',3);
        %nMeno_prefHead_smooth = smoothdata(nMeno_binAngle,'gaussian',3);
        
%         figure();
%         plot(centers_angle,nMeno_binAngle);
%         hold on
%         plot(centers_angle,nMeno_prefHead_smooth)
        
%         if sum(~isnan(Meno_binAngle)) < length(Meno_binAngle)  - 2 
%             Meno_prefHead = nan; 
%             Meno_prefHead_amp = nan;
%             
%         else
%             
%             [Meno_bump_params, model_data_meno, x_grid_meno] = fit_sinusoid(Meno_prefHead_smooth, [0,2*pi], 0);
%             Meno_prefHead = wrapTo180((rad2deg(wrapToPi(Meno_bump_params.pos_rad + deg2rad(15)))) - 180);
%             
%             %if Meno_bump_params.adj_rs > 0.4
%                 Meno_prefHead = Meno_prefHead;
                Meno_prefHead_amp = max(Meno_prefHead) - min(Meno_prefHead); %Meno_bump_params.amp;
%             else
%                 Meno_prefHead = nan;
%                 Meno_prefHead_amp = nan;
%             end
            
%         end
        
%         if sum(~isnan(nMeno_binAngle)) < length(nMeno_binAngle)  - 2 
%             nMeno_prefHead = nan;
%             nMeno_prefHead_amp = nan;
%         else
%             [nMeno_bump_params, model_data_nMeno, x_grid_nMeno] = fit_sinusoid(nMeno_prefHead_smooth, [0,2*pi], 0);
%             nMeno_prefHead = wrapTo180((rad2deg(wrapToPi(nMeno_bump_params.pos_rad + deg2rad(15)))) - 180);
%             
%             %if nMeno_bump_params.adj_rs > 0.5
%                 nMeno_prefHead = nMeno_prefHead;
                nMeno_prefHead_amp = max(nMeno_prefHead) - min(nMeno_prefHead); %nMeno_bump_params.amp;
%             else
%                 nMeno_prefHead = nan;
%                 nMeno_prefHead_amp = nan;
%             end
%         end
        
%         bump_params = fit_sinusoid(prefHead_smooth, [0,2*pi], 1);
%         prefHead = wrapTo180((rad2deg(wrapToPi(bump_params.pos_rad + deg2rad(15)))) - 180);
%         figure();clf;
%         plot(centers_angle,sum_mean{2})
%         hold on

% allGoals = repelem(MenoDataFly.Goal,round(MenoDataFly.timeMov)); 
% if size(allGoals,2) > size(allGoals,1)
%     allGoals = allGoals';
% end
% ave_goal = -(circ_mean(allGoals)); % switch from cue pos to HD
% avePrefHead = -mean(prefHead,'omitnan');
% 
% goalDiff = rad2deg(angdiff(ave_goal,deg2rad(avePrefHead))); 
% 
% figure(11);clf;
% set(gcf,'color','w','renderer','painters')
% a1 = subplot(1,1,1);
% hold on
% plot(-centers_angle,Meno_prefHead,'k')
% plot(-centers_angle,nMeno_prefHead,'color',[0.75,0.75,0.75])
% xlabel('HD')
% ylabel('Vm')
% 
% if regexp(rootDir,'PFL3')
%     title([cellSide,' goalDiff: ',num2str(goalDiff)])
% else
%     title('goalDiff: ',num2str(goalDiff))
% end
% box off
% xlim(a1,[-180,180])
% maxMeno = max(Meno_prefHead);
% minMeno = min(Meno_prefHead);

% fig2 = figure();
% set(gcf,'color','w','renderer','painters')
% b1 = subplot(1,1,1);
% plot(-centers_angle,nMeno_prefHead,'k')
% xlabel('HD')
% ylabel('Vm')
% title('nMeno')
% box off
% maxnMeno = max(nMeno_prefHead);
% minnMeno = min(nMeno_prefHead);

%ylim(a1,[min(minMeno,minnMeno),max(maxMeno,maxnMeno)])

%ylim(b1,[min(minMeno,minnMeno),max(maxMeno,maxnMeno)])
%xlim(b1,[-180,180])
%set(fig2,'position',[-800,600,500,250])
%set(fig1,'position',[-1500,600,500,250])
% 
% saveas(fig1,fullfile('Z:\Dropbox (HMS)\Wilson_Lab_Data\Analysis\Westeinde2022_ms\EphysPlots\preferredHeadings\meno',['menoPrefHead_fly',num2str(i),'.fig']))
% saveas(fig1,fullfile('Z:\Dropbox (HMS)\Wilson_Lab_Data\Analysis\Westeinde2022_ms\EphysPlots\preferredHeadings\meno',['menoPrefHead_fly',num2str(i),'.svg']))
% 
% saveas(fig2,fullfile('Z:\Dropbox (HMS)\Wilson_Lab_Data\Analysis\Westeinde2022_ms\EphysPlots\preferredHeadings\nMeno',['nMenoPrefHead_fly',num2str(i),'.fig']))
% saveas(fig2,fullfile('Z:\Dropbox (HMS)\Wilson_Lab_Data\Analysis\Westeinde2022_ms\EphysPlots\preferredHeadings\nMeno',['nMenoPrefHead_fly',num2str(i),'.svg']))
%         
        
        
        
        
        
        Meno_Sum(i) = Meno_prefHead_amp;
        nMeno_Sum(i) = nMeno_prefHead_amp;

%         summaryArray(count).fly = fly;
%         summaryArray(count).Meno_prefHead_amp =  Meno_prefHead_amp;
%         summaryArray(count).Meno_prefHead_amp =  Meno_prefHead;
%         summaryArray(count).nMeno_prefHead_amp =  nMeno_prefHead_amp;
%         summaryArray(count).nMeno_prefHead_amp =  nMeno_prefHead;

end  
% lineplots    

                
                figure()
                set(gcf,'color','w')
                set(gcf,'renderer','painters')
                hold on
                allSum = [nMeno_Sum',Meno_Sum'];
                for fly = 1:size(allSum,1)                
                    plot([1,1.5],allSum(fly,:),'-o','color',[0.75,0.75,0.75],'MarkerFaceColor',[0.75,0.75,0.75]);
                end
                plot([1,1.5],mean(allSum,1),'k-o','MarkerFaceColor','k')
                xlim([0.5,2])
                [pVm_wil, hVm_wil] = signrank(allSum(:,1),allSum(:,2));
                [htt, ptt] = ttest(allSum(:,2),allSum(:,1));
                ylabel('prefHead amp')
                xlabel('not meno vs meno')
                

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

    if savePlots == 1
        saveas(z, fullfile(lineplotDir,[expID,'_',num2str(nTrial),'_zScore_behaviour_no0vel.fig']));
    end
end