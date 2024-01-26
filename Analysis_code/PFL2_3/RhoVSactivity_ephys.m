    folders_PFL2 = get_folders_ephys('Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL2\new_flies');
    folders_PFL3 = get_folders_ephys('Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL3\new_flies');

    all_data = [folders_PFL2;folders_PFL3];
     
     all_folders = []; 
     for f = 1:length(all_data)
        all_folders{f} = (all_data(f).folder);
     end
     all_folders = all_folders';
  folders = [];    
   cnt = 0;
   %all_folders = summaryArray.Folder;
    for f = 1:size(all_folders,1)
       %if isempty(regexp(all_folders{f},'R'))
            cnt = cnt + 1;
            folders{cnt} = all_folders{f};
            slash = regexp(all_folders{f},'\');
            expIDidx = [slash(6)+1:slash(7)-1];
            expID = char(all_folders{f});
            expID = expID(expIDidx);
            flies(cnt) = string(expID);
%         [start,finish] = regexp(all_folders(f),'flies\', 'ignorecase');
%         dateIdx = regexp(all_folders(f),'\');
%         dateIdx = [dateIdx(3)+1:dateIdx(4)-1]; 
%         Date = char(all_folders(f));
%         Date = Date(dateIdx);
%         if isempty(finish)
%             [start,finish] = regexp(all_folders(f),'_Fly');
%         end
%         fly_temp = char(all_folders(f));
%         fly = fly_temp(start:finish + 1);
%         flyID = strcat(Date,fly);
%         flies(f) = string(flyID);
       %end
    end
        flies = flies';
        uniqueFlies = unique(flies); 
        flyCount = zeros(size(folders));
        for fly = 1:length(uniqueFlies)
            flyCount(flies == uniqueFlies(fly)) = fly;
        end
        
        %
flySum = []; 
flySum_var = []; 
flySum_3Dnum = []; 
flySum_3D = []; 
flySum_angle = [];
flySum__angleVar = [];
for currentFly = 1:length(uniqueFlies) % goes through each fly
    flyTrials = folders(flyCount == currentFly)';
    for trial = 1:size(flyTrials,1) % goes through each trial for a fly
        %try
       % folder = table2array(summaryArray(trial,1)); 
       folder = flyTrials{trial};
        if strcmp(folder(end),'.')
            folder = folder(1:end-2); 
        end
        
        processedDir = fullfile(folder,'processedData');
        clear bData tData 
        try
            load(fullfile(processedDir,'pro_behaviourData.mat'))
            bData = pro_behaviourData{1};
            load(fullfile(processedDir,'pro_trialData.mat'))
            tData = pro_trialData{1};
        catch
            load(fullfile(folder,'pro_behaviourData.mat'))
            bData = processed_behaviourData{1};
            load(fullfile(folder,'pro_trialData.mat'))
            tData = processed_trialData{1};
        end

        velFor = [];
        velYaw = []; 
        velSide = []; 
        VM = []; 
        FR = []; 
        rho = [];
        theta = []; 
        cueAngle = [];
        
        window = 30; 
        step = 0.1 * 1000; 
        minVel = 3; 

        window = window * 1000; 
        if window > length(bData.vel_for)
            window = length(bData.vel_for); 
        end
        count = 1;
        jump_idx = find(bData.jumps == 1);
        
        no0vel_idx = find(abs(bData.vel_for) + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw)*4.5) > minVel);
        angle = bData.angle(no0vel_idx);
        all_vf = bData.vel_for(no0vel_idx);
        all_vy = bData.vel_yaw(no0vel_idx);
        all_vs = bData.vel_side(no0vel_idx);
        all_FR = tData.fRate_sec(no0vel_idx);
        try
            all_vm = tData.smoothVm(no0vel_idx);
        catch
            all_vm = tData.smooth_Vm(no0vel_idx);
        end


        for i = step/2:step:length(angle) - step/2

            idx = i - window/2:i + window/2; 

            if idx(end) > length(angle)-step/2 && idx(1) < step/2 
                idx = step/2:1:length(angle)-step/2;
            elseif idx(end) > length(angle)-step/2
                idx = idx(1):1:length(angle)-step/2;
            elseif idx(1) < step/2 
                idx = step/2:idx(end); 
            end            

            if sum(ismember(idx, jump_idx)) % remove idx from window that were influenced by jumps
                badIdx = ismember(idx, jump_idx);
                idx = idx(badIdx == 0);
            end

                angles_flyFor = angle(idx); 
                if ~isempty(angles_flyFor)
                    x = cosd(angles_flyFor); 
                    y = sind(angles_flyFor); %my arena has - angles to the left of the fly, + to the right, multiply y component by -1 to align physical arena coords to polar plot angles
                    idx_windows{count,1} = idx;
                    mean_headingVectors(1)= sum(x)/length(x); 
                    mean_headingVectors(2)= sum(y)/length(y);
                else
                    mean_headingVectors(1)= nan; 
                    mean_headingVectors(2)= nan; 
                    idx_windows{count,1} = idx;
                    slope(count) = nan;
                    count = count + 1;
                end

            velFor(count) = all_vf(i); 
            velYaw(count) = all_vy(i);
            velSide(count) = all_vs(i);
            VM(count) = all_FR(i); 
            FR(count) = all_FR(i);
            cueAngle(count) = angle(i);
            rho(count) = sqrt(mean_headingVectors(1)^2 + mean_headingVectors(2)^2);
            theta(count) = atan2(mean_headingVectors(2),mean_headingVectors(1)); 
            
%             
%             if rho(count) > 0.9
%                 figure();
%                 plot(angle(idx))
%                 pause
%             end
            count = count + 1;

        end
        
        edges_rho = [0:0.1:1];
        [VMBinned, ~, variance] = binData(VM', rho', edges_rho);
        
        flySum(currentFly,trial,:) = VMBinned;
        flySum_var(currentFly,trial,:) = variance;
        
        [angleBinned, centers_rho, variance] = binData(cueAngle', rho', edges_rho);
        
        flySum_angle(currentFly,trial,:) = angleBinned;
        flySum_angleVar(currentFly,trial,:) = variance;
        
        edges_angle = [-180:20:180];
        
        [N, rho_cueAngle_vm, x_centers, y_centers] = create_activity_heatmap(rho', cueAngle', VM', edges_rho, edges_angle);
        
        flySum_3D(currentFly,trial,:,:) = rho_cueAngle_vm;
        flySum_3Dnum(currentFly,trial,:,:) = N;

           
    end
end

%% angle
nCol = length(uniqueFlies);
colormap = cbrewer2('RdBu',nCol);

count = 1; 
nonZeroTrials = [];
 figure();hold on;
for currentFly = 1:length(uniqueFlies)
   trials = squeeze(flySum_angle(currentFly,:,:));
    for t = 1:size(trials,1)
        trialVM = trials(t,:); 
        if sum(trialVM,'omitnan') ~= 0 
            %trialVM = trialVM - trialVM(1);
            nonZeroTrials(count,:) = trialVM; 
            count = count + 1; 
            scatter(centers_rho, trialVM,[],colormap(currentFly,:))
        end
    end        
end
 scatter(centers_rho,mean(nonZeroTrials,1,'omitnan'),'k','filled')


nonZeroTrials = [];
figure();hold on;
count = 1; 
for currentFly = 1:length(uniqueFlies)
   trials = squeeze(flySum_angleVar(currentFly,:,:));
    for t = 1:size(trials,1)
        trialVM = trials(t,:); 
        if sum(trialVM,'omitnan') ~= 0 
            nonZeroTrials(count,:) = trialVM; 
            %trialVM = trialVM - trialVM(1);
            scatter(centers_rho, nonZeroTrials(count,:),[],colormap(currentFly,:))
            count = count + 1; 
        end
    end        
end
 scatter(centers_rho,mean(nonZeroTrials,1,'omitnan'),'k','filled')
 
 %% activity
 
 nCol = length(uniqueFlies);
colormap = cbrewer2('RdBu',nCol);

count = 1; 
nonZeroTrials = [];
figure();hold on;
for currentFly = 1:length(uniqueFlies)
   trials = squeeze(flySum(currentFly,:,:));
    for t = 1:size(trials,1)
        trialVM = trials(t,:); 
        if sum(trialVM,'omitnan') ~= 0 
            %trialVM = trialVM - trialVM(1);
            nonZeroTrials(count,:) = trialVM; 
            count = count + 1; 
            scatter(centers_rho, trialVM,[],colormap(currentFly,:))
        end
    end        
end
 scatter(centers_rho,mean(nonZeroTrials,1,'omitnan'),'k','filled')
xlabel('rho')
ylabel('Vm')

nonZeroTrials = [];
figure();hold on;
count = 1; 
for currentFly = 1:length(uniqueFlies)
   trials = squeeze(flySum_var(currentFly,:,:));
    for t = 1:size(trials,1)
        trialVM = trials(t,:); 
        if sum(trialVM,'omitnan') ~= 0 
            nonZeroTrials(count,:) = trialVM; 
            count = count + 1; 
            %trialVM = trialVM - trialVM(1);
            scatter(centers_rho, trialVM,[],colormap(currentFly,:))
        end
    end        
end
 scatter(centers_rho,mean(nonZeroTrials,1,'omitnan'),'k','filled')
 xlabel('rho')
ylabel('Vm variance')
 
 %% heatmaps
 
 for currentFly = 1:length(uniqueFlies)
   trials = squeeze(flySum_3D(currentFly,:,:,:));
    for t = 1:size(trials,1)
        temp = squeeze(flySum_3D(currentFly,t,:,:));
        temp_count = squeeze(flySum_3Dnum(currentFly,t,:,:));
        temp(temp_count == 0) = nan;
        flySum_3D(currentFly,t,:,:) = temp;        
        trial_HD_var(currentFly,t,:) = var(temp,0,1,'omitnan');         
    end
 end

 

figure();hold on;
all_nonZeroTrials = [];
count2 = 1; 
for currentFly = 1:length(uniqueFlies)
   nonZeroTrials = [];
   count = 1; 
   trials = squeeze(trial_HD_var(currentFly,:,:));
    for t = 1:size(trials,1)
        trialVM = trials(t,:); 
        if sum(trialVM,'omitnan') ~= 0 
            nonZeroTrials(count,:) = trialVM; 
            all_nonZeroTrials(count2,:) = trialVM; 
            count2 = count2 + 1; 
            count = count + 1; 
            %trialVM = trialVM - trialVM(1);
           scatter(centers_rho, trialVM,[],colormap(currentFly,:))
        end
    end 
    %scatter(centers_rho,mean(nonZeroTrials,1,'omitnan'),[],colormap(currentFly,:))
end
scatter(centers_rho,mean(all_nonZeroTrials,1,'omitnan'),'k','filled')
xlabel('rho')
ylabel('Vm variance across HD')
 
 
aveTrial_3D = squeeze(mean(flySum_3D,2,'omitnan')); 
aveFly_3D = squeeze(mean(aveTrial_3D,1,'omitnan')); 

for fly = 1:size(aveTrial_3D,1)
    trial_map = squeeze(aveTrial_3D(fly,:,:));
    ncol = length(unique(trial_map));
    color = flipud(cbrewer2('RdYlBu', ncol));
    
    figure();
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    s = pcolor(trial_map);
    colormap(color)
    xt = linspace(1,numel(x_centers),5);                            
    xtlbl = linspace(x_centers(xt(1)), x_centers(xt(end)), 5);   
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
    yt = linspace(1,numel(y_centers),5); 
    ytlbl = linspace(y_centers(yt(1)),y_centers(yt(end)), 5);
    set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
    set(s, 'EdgeColor', 'none');
    set(gca,'color','none')
    colorbar
    box off
end
