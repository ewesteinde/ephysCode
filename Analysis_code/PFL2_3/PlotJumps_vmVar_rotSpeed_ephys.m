function PlotJumps_vmVar_rotSpeed_ephys(rootDir)       

     folders1 = get_folders_ephys(rootDir);
     folders2 = get_folders_ephys('Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL2\new_flies');
     
     folders = [folders1; folders2];
     
     headers = {'folder','trial','jumpSize','jumpIdx','velFor','velSide','velYaw','cueAngle','corrected','timeMov','prefHead','preJump','postJump','jump_vm','cueDiff'};
     jump_summary = cell2table(cell(1,15),'VariableNames',headers); 
     count = 1; 

    %% Process each folder
    folderNum = length(folders);
    fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
    countFail = 1;
    for ff = 1:folderNum
      
      %% Get folder information
      folder = folders(ff).folder;
      
       
%% Load in fictrac & ROI data
       %try
            if strcmp(folder(end),'.')
                folder = folder(1:end-2); 
            end
            
            try
                processedDir = fullfile(folder,'processedData');
                load(fullfile(processedDir,'pro_behaviourData.mat'))
                load(fullfile(processedDir,'pro_trialData.mat'))
                bData = pro_behaviourData{1}; 
                tData = pro_trialData{1}; 
            catch
                load(fullfile(folder,'pro_behaviourData.mat'))
                load(fullfile(folder,'pro_trialData.mat'))
                bData = processed_behaviourData{1}; 
                tData = processed_trialData{1};  
            end

        numTrials = 1;

        %% visualize individual trials
        jumpWin = 15;

        for t = 1:numTrials

        [jump_array_rho, ~, ~] = detect_jumps_ephys(bData.frY, 30,jumpWin, 1000);
        [jump_array, ~, ~] = detect_jumps_ephys(bData.frY, jumpWin,jumpWin, 1000);

        total_mov_mm = abs(bData.vel_for + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw)*4.5));
        
        try
            %activity = (tData.smooth_Vm - median(tData.smooth_Vm))/mad(tData.smooth_Vm);
            activity = tData.smooth_Vm;
        catch
            %activity = (tData.smoothVm - median(tData.smoothVm))/mad(tData.smoothVm);
            activity = tData.smoothVm;
        end
        
%        activity = tData.fRate_sec;

        count90 = 1;
        count180 = 1;
        jump90_time = []; 
        jump180_time = [];
        
        jumpIdx = [jump_array(1,1):jump_array(1,3)]';  
        jumpIdx_rho = [jump_array_rho(1,1):jump_array_rho(1,3)]';
        if jumpIdx(1) < 1 || jumpIdx_rho(1) < 1
            jump_array(1,:) = []; 
            jump_array_rho(1,:) = []; 
        end
        
        for jump = 1:size(jump_array,1)
            jumpIdx = [jump_array(jump,1):jump_array(jump,3)]'; 
            if jumpIdx(end) > size(bData,1)
                jump_array(jump,:) = []; 
                jump_array_rho(jump,:) = []; 
            end
        end
                
                for jump = 1:size(jump_array,1)
                    jumpSize = jump_array(jump,4);
                    jumpIdx = [jump_array(jump,1):jump_array(jump,3)]'; 
                    jump_vf = bData.vel_for(jumpIdx);
                    jump_vy = bData.vel_yaw(jumpIdx);
                    jump_vs = bData.vel_side(jumpIdx);
                    jump_angle = bData.angle(jumpIdx);

                    preJumpIdx = [jump_array(jump,1):jump_array(jump,2)-1];
                    preJumpIdx_rho = [jump_array_rho(jump,1):jump_array_rho(jump,2)-1];
                    timeMov = (sum(abs(bData.vel_for(preJumpIdx_rho)) + abs(bData.vel_side(preJumpIdx_rho)) + abs(deg2rad(bData.vel_yaw(preJumpIdx_rho))*4.5) > 3))/1000; % seconds
                    preJump_base = mean(activity(preJumpIdx));
                    jump_vm = activity(jumpIdx) - preJump_base; 
                    

                    [rho, ~] = CalculateAverageHeading_ephys(bData,1.5, preJumpIdx_rho);
                    [~, preJumpAngle] = CalculateAverageHeading_ephys(bData,0, preJumpIdx);
                    
                    %preJumpAngle = ftT.cueAngle{1}(preJumpIdx(end));
                    
                    postJumpIdx = jump_array(jump,2)+1;
                    sizePostJump = size(bData.angle(postJumpIdx:jump_array(jump,3)));
                    preJump = ones(sizePostJump) * preJumpAngle;
%                     preJump_head = ones(sizePostJump) * prefHead;
%                     
%                     %diff b/w pref head & post jump angle
%                     cueDiff_head = angdiff(preJump_head,deg2rad(bData.angle(postJumpIdx:jump_array(jump,3)))); 
%                     to_prefHead = sum(abs(cueDiff_head) < deg2rad(30));
%                     away = 0; 
%                     towards = 0; 
%                     
%                     if abs(rad2deg(angdiff(preJumpAngle,prefHead))) < 30
%                         away = 1; 
%                     elseif to_prefHead > 0     
%                         towards = 1; 
%                     end
                    
                    cueDiff = angdiff(preJump,deg2rad(bData.angle(postJumpIdx:jump_array(jump,3)))); 
                    
                    if regexp(folder,'070121')
                        trueJumpSize = rad2deg(angdiff(deg2rad(bData.angle(jump_array(jump,2) - 100)),deg2rad(bData.angle(jump_array(jump,2) + 100)))); 
                        t1 = 90-trueJumpSize;
                        t2 = -90 - trueJumpSize;
                        
                        if t1 < t2
                            jumpSize = 90;
                        else
                            jumpSize = -90; 
                        end
                    
                    end
                    
                    if jumpSize == 180
                        correct = sum(abs(cueDiff) < deg2rad(45),'omitnan'); 
                    else
                        correct = sum(abs(cueDiff) < deg2rad(30),'omitnan');
                    end 
                        
                    
                    jump_summary.folder(count) = {folder}; 
                    jump_summary.trial(count) = {t}; 
                    jump_summary.jumpSize(count) = {jumpSize}; 
                    jump_summary.jumpIdx(count) = {jumpIdx};
                    jump_summary.velFor(count) = {jump_vf};
                    jump_summary.velYaw(count) = {jump_vy};
                    jump_summary.velSide(count) = {jump_vs}; 
                    jump_summary.cueAngle(count) = {jump_angle};
                    jump_summary.timeMov(count) = {timeMov}; 
                    jump_summary.preJump(count) = {bData.angle(jump_array(jump,2)-1)}; 
                    jump_summary.postJump(count) = {bData.angle(jump_array(jump,2)+1)}; 
                    jump_summary.prefHead(count) = {[]};%{prefHead}; 
                    jump_summary.jump_vm(count) = {jump_vm};
                    jump_summary.cueDiff(count) = {cueDiff}; 

                  
                    if  timeMov >= 15  && correct >= 10 && rho >= 0.7
                        jump_summary.corrected(count) = {1};
                        %figure(1);clf;plot(jump_angle)
                    elseif correct >= 10 && rho < 0.7
                        jump_summary.corrected(count) = {2}; 
                        %figure(1);clf;plot(jump_angle)
                    elseif correct <= 0 && timeMov >= 15
                        jump_summary.corrected(count) = {0}; 
                        %figure(1);clf;plot(jump_angle)
                    else
                        jump_summary.corrected(count) = {3}; 
                        %figure(1);clf;plot(jump_angle)
                        % corrected, rho > threshold, timeMov < threshold
                    end
                    
                    count = count + 1; 
                    
                end
        end
%         for jump = 1:size(jump_summary,1)
%             figure();
%             plot(jump_summary.cueAngle{jump})
%         end

%         catch
%             disp([folder, ' failed'])
%         end
    end  
    Sum180 = jump_summary(cell2mat(jump_summary.jumpSize) == 180 & cell2mat(jump_summary.corrected) == 1,:); 
    nSum180 = jump_summary(cell2mat(jump_summary.jumpSize) == 180 & cell2mat(jump_summary.corrected) == 0,:); 
    Sum90 = jump_summary(cell2mat(jump_summary.jumpSize) == 90 & cell2mat(jump_summary.corrected) == 1,:);
    nSum90 = jump_summary(cell2mat(jump_summary.jumpSize) == 90 & cell2mat(jump_summary.corrected) == 0,:);
    Summ90 = jump_summary(cell2mat(jump_summary.jumpSize) == -90 & cell2mat(jump_summary.corrected) == 1,:);
    nSumm90 = jump_summary(cell2mat(jump_summary.jumpSize) == -90 & cell2mat(jump_summary.corrected) == 0,:);
    %%
    jumpTime = [-jumpWin+(1/1000):1/1000:jumpWin+(1/1000)];
    plotWin = 2; 
    vfSum = zeros(size(Sum180,1),size(jump_vf,1));
    vySum = zeros(size(Sum180,1),size(jump_vf,1));
    vsSum = zeros(size(Sum180,1),size(jump_vf,1));
    VmSum = zeros(size(Sum180,1),size(jump_vf,1));
    cueDiff_sum = zeros(size(Sum180,1),size(cueDiff,1));
    figure();
    sgtitle('180 jumps')
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    ax1 = subplot(3,1,1);
    hold on
    ax2 = subplot(3,1,2);
    hold on
    ax3 = subplot(3,1,3); 
    hold on
%     ax4 = subplot(3,1,3); 
%     hold on
    for jump = 1:size(Sum180,1)  
        %plot(ax1, jumpTime,Sum180.cueAngle{jump})
        %plot(ax2,jumpTime,Sum180.velFor{jump})
        VmSum(jump,:) = Sum180.jump_vm{jump};
        vfSum(jump,:) = Sum180.velFor{jump}'; 
        %plot(ax3,jumpTime,abs(Sum180.velYaw{jump}))
        vySum(jump,:) = Sum180.velYaw{jump}';
        %plot(ax4,jumpTime,abs(Sum180.velSide{jump}))
        vsSum(jump,:) = Sum180.velSide{jump}';
        cueDiff_sum(jump,:) = Sum180.cueDiff{jump}';
    end
    
    
    SEM_Vm = std(VmSum,[],1,'omitnan') / sqrt(size(VmSum,1));
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vs = std(vsSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    nvfSum = zeros(size(nSum180,1),size(jump_vf,1));
    nvySum = zeros(size(nSum180,1),size(jump_vf,1));
    nvsSum = zeros(size(nSum180,1),size(jump_vf,1));
    nVmSum = zeros(size(nSum180,1),size(jump_vf,1));
    ncueDiff_sum = zeros(size(nSum180,1),size(cueDiff,1));
    
    for jump = 1:size(nSum180,1)  
        %plot(ax1, jumpTime,Sum180.cueAngle{jump})
        %plot(ax2,jumpTime,Sum180.velFor{jump})
        nVmSum(jump,:) = nSum180.jump_vm{jump};
        nvfSum(jump,:) = nSum180.velFor{jump}'; 
        %plot(ax3,jumpTime,abs(Sum180.velYaw{jump}))
        nvySum(jump,:) = nSum180.velYaw{jump}';
        %plot(ax4,jumpTime,abs(Sum180.velSide{jump}))
        nvsSum(jump,:) = nSum180.velSide{jump}';
        ncueDiff_sum(jump,:) = nSum180.cueDiff{jump}';
    end
    nSEM_Vm = std(nVmSum,[],1,'omitnan') / sqrt(size(nVmSum,1));
    nSEM_vf = std(nvfSum,[],1,'omitnan') / sqrt(size(nvfSum,1));
    nSEM_vy = std(nvySum,[],1,'omitnan') / sqrt(size(nvfSum,1));
    nSEM_vs = std(nvsSum,[],1,'omitnan') / sqrt(size(nvfSum,1));
    
% 
%     figure();
%     sgtitle('180 jumps')
%     set(gcf,'color','w')
%     set(gcf,'renderer','painters')
%     ax1 = subplot(1,1,1);
%     hold on
    keepIndex = ~isnan(nSEM_Vm);
    SEMhigh = [smoothdata(mean(abs(nVmSum(:,keepIndex)),'omitnan'),'gaussian',1) + nSEM_Vm(keepIndex)]; 
    SEMlow = [smoothdata(mean(abs(nVmSum(:,keepIndex)),'omitnan'),'gaussian',1) - nSEM_Vm(keepIndex)];
    %SEMhigh = [smoothdata(var(VmSum(:,keepIndex),'omitnan'),'gaussian',1) + SEM_Vm(keepIndex)]; 
    %SEMlow = [smoothdata(var(VmSum(:,keepIndex),'omitnan'),'gaussian',1) - SEM_Vm(keepIndex)];
    %plot(ax1,jumpTime,VmSum,'color',[0.5,0.5,0.5],'LineWidth',0.5)
     plot(ax1,jumpTime,(mean(abs(nVmSum),1)),'r','LineWidth',1.5)
    %plot(ax1,jumpTime,((nVmSum)),'color',[0,0,0,0.2])
   % plot(ax1,jumpTime,var(VmSum,1),'k','LineWidth',1.5)
    patch(ax1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0,0],'FaceAlpha',0.25,'EdgeColor','none')
    xline(ax1,0)
    xlim(ax1,[-plotWin,plotWin])
    ylabel(ax1,'|Z(Vm)|')
    xlabel(ax1,'Time (s)')
    
%     figure();
%     sgtitle('180 jumps')
%     set(gcf,'color','w')
%     set(gcf,'renderer','painters')
%     ax1 = subplot(1,1,1);
%     hold on
    keepIndex = ~isnan(SEM_Vm);
    SEMhigh = [smoothdata(mean(abs(VmSum(:,keepIndex)),'omitnan'),'gaussian',1) + SEM_Vm(keepIndex)]; 
    SEMlow = [smoothdata(mean(abs(VmSum(:,keepIndex)),'omitnan'),'gaussian',1) - SEM_Vm(keepIndex)];
    %SEMhigh = [smoothdata(var(VmSum(:,keepIndex),'omitnan'),'gaussian',1) + SEM_Vm(keepIndex)]; 
    %SEMlow = [smoothdata(var(VmSum(:,keepIndex),'omitnan'),'gaussian',1) - SEM_Vm(keepIndex)];
    %plot(ax1,jumpTime,VmSum,'color',[0.5,0.5,0.5],'LineWidth',0.5)
    plot(ax1,jumpTime,(mean(abs(VmSum),1)),'k','LineWidth',1.5)
     %plot(ax1,jumpTime,((VmSum)),'color',[0,0,0,0.1])
    %plot(ax1,jumpTime,var(VmSum,1),'k','LineWidth',1.5)
    patch(ax1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xline(ax1,0)
    xlim(ax1,[-plotWin,plotWin])
    ylabel(ax1,'|Z(Vm)|')
    xlabel(ax1,'Time (s)')
    
    %xlim(ax1,[min(jumpTime),max(jumpTime)])
    keepIndex = ~isnan(SEM_vf);
    SEMhigh = [mean(vfSum(:,keepIndex),'omitnan') + SEM_vf(keepIndex)]; 
    SEMlow = [mean(vfSum(:,keepIndex),'omitnan') - SEM_vf(keepIndex)];
    %plot(ax2,jumpTime,vfSum,'color',[0.5,0.5,0.5],'LineWidth',0.5)
    plot(ax3,jumpTime,mean(vfSum,1),'k','LineWidth',1.5)
    patch(ax3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xline(ax3,0)
    xlim(ax3,[-plotWin,plotWin])
    ylabel(ax3,'forward velocity (mm/s)')
    xlabel(ax3,'Time (s)')
    
    %xlim(ax1,[min(jumpTime),max(jumpTime)])
    keepIndex = ~isnan(nSEM_vf);
    SEMhigh = [mean(nvfSum(:,keepIndex),'omitnan') + nSEM_vf(keepIndex)]; 
    SEMlow = [mean(nvfSum(:,keepIndex),'omitnan') - nSEM_vf(keepIndex)];
    %plot(ax2,jumpTime,vfSum,'color',[0.5,0.5,0.5],'LineWidth',0.5)
    plot(ax3,jumpTime,mean(nvfSum,1),'r','LineWidth',1.5)
    patch(ax3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0,0],'FaceAlpha',0.25,'EdgeColor','none')
    xline(ax3,0)
    xlim(ax3,[-plotWin,plotWin])
    ylabel(ax3,'forward velocity (mm/s)')
    xlabel(ax3,'Time (s)')
    
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [mean(abs(vySum(:,keepIndex)),'omitnan') + SEM_vy(keepIndex)]; 
    SEMlow = [mean(abs(vySum(:,keepIndex)),'omitnan') - SEM_vy(keepIndex)];
    %plot(ax3,jumpTime,abs(vySum),'color',[0.5,0.5,0.5],'LineWidth',0.5)
    plot(ax2,jumpTime,mean(abs(vySum),1),'k','LineWidth',1.5)
    patch(ax2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xline(ax2,0)
    xlim(ax2,[-plotWin,plotWin])
    ylabel(ax2,'rotational speed (deg/s)')
    xlabel(ax2,'Time (s)')
    
    keepIndex = ~isnan(nSEM_vy);
    SEMhigh = [mean(abs(nvySum(:,keepIndex)),'omitnan') + nSEM_vy(keepIndex)]; 
    SEMlow = [mean(abs(nvySum(:,keepIndex)),'omitnan') - nSEM_vy(keepIndex)];
    %plot(ax3,jumpTime,abs(vySum),'color',[0.5,0.5,0.5],'LineWidth',0.5)
    plot(ax2,jumpTime,mean(abs(nvySum),1),'r','LineWidth',1.5)
    patch(ax2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0,0],'FaceAlpha',0.25,'EdgeColor','none')
    xline(ax2,0)
    xlim(ax2,[-plotWin,plotWin])
    ylabel(ax2,'rotational speed (deg/s)')
    xlabel(ax2,'Time (s)')

linkaxes([ax1,ax2,ax3],'x')


figure; 
sgtitle('180 jumps')
set(gcf,'color','w')
set(gcf,'renderer','painters')
yyaxis left
hold on
keepIndex = ~isnan(SEM_Vm);
SEMhigh = [mean(abs(VmSum(:,keepIndex)),'omitnan') + SEM_Vm(keepIndex)]; 
SEMlow = [mean(abs(VmSum(:,keepIndex)),'omitnan') - SEM_Vm(keepIndex)];
patch([jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,1],'FaceAlpha',0.1,'EdgeColor','none')
plot(jumpTime,(mean(abs(VmSum),1)),'b','LineWidth',1.5)
xline(0)
%ylim([min(SEMlow),max(SEMhigh)])
ylim([0.5,8])
ylabel('|Z(Vm)|')
yyaxis right
hold on
keepIndex = ~isnan(SEM_vy);
SEMhigh = [mean(abs(vySum(:,keepIndex)),'omitnan') + SEM_vy(keepIndex)]; 
SEMlow = [mean(abs(vySum(:,keepIndex)),'omitnan') - SEM_vy(keepIndex)];
patch([jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[1,0,0],'FaceAlpha',0.1,'EdgeColor','none')
plot(jumpTime,mean(abs(vySum),1),'r','LineWidth',1.5)
%ylim([min(SEMlow),max(SEMhigh)])
ylim([
ylabel('|Yaw|')
xlim([-0.2,1])

% figure(); 
% hold on
% for j = 1:size(jump_summary(cell2mat(jump_summary.jumpSize) == 180,:).cueDiff,1)
%     a = plot(jumpTime(30002:end), rad2deg(jump_summary(cell2mat(jump_summary.jumpSize) == 180,:).cueDiff{j}), 'color',[0.25,0.25,0.25],'LineWidth',0.5);
%     a.YData(abs(diff(a.YData)) > pi) = nan;
% end
% yline(-60,'r')
% yline(60,'r')
% yline(-30,'b')
% yline(30,'b')
% title('180')
% 
% figure(); 
% hold on
% for j = 1:size(jump_summary(abs(cell2mat(jump_summary.jumpSize)) == 90,:).cueDiff,1)
%     a = plot(jumpTime(30002:end), rad2deg(jump_summary(abs(cell2mat(jump_summary.jumpSize)) == 90,:).cueDiff{j}), 'color',[0.25,0.25,0.25],'LineWidth',0.5);
%     a.YData(abs(diff(a.YData)) > pi) = nan;
% end
% yline(-60,'r')
% yline(60,'r')
% yline(-30,'b')
% yline(30,'b')
% title('90')

 figure;
 set(gcf,'renderer','painters','color','w')
subplot(2,1,1);
sgtitle('180 jumps')
plot(jumpTime,VmSum,'k')
xlim([-plotWin,plotWin])
box off
subplot(2,1,2)
plot(jumpTime,nVmSum,'r')
xlim([-plotWin,plotWin])
box off

        
    vfSum = zeros(size(Sum90,1),size(jump_vf,1));
    vySum = zeros(size(Sum90,1),size(jump_vf,1));
    vsSum = zeros(size(Sum90,1),size(jump_vf,1));
    VmSum = zeros(size(Sum90,1),size(jump_vf,1));
    cueDiff_sum = zeros(size(Sum90,1),size(cueDiff,1));
    figure();
    sgtitle('90 jumps')
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    bx1 = subplot(3,1,1); 
    hold on
    bx2 = subplot(3,1,2);
    hold on
    bx3 = subplot(3,1,3);
    hold on
%     bx4 = subplot(3,1,3); 
%     hold on
    for jump = 1:size(Sum90,1)  
        %plot(bx1, jumpTime,Sum90.cueAngle{jump})
        %plot(bx2,jumpTime,Sum90.velFor{jump})
        VmSum(jump,:) = Sum90.jump_vm{jump};
        vfSum(jump,:) = Sum90.velFor{jump}'; 
        %plot(bx3,jumpTime,Sum90.velYaw{jump})
        vySum(jump,:) = Sum90.velYaw{jump}';
        %plot(bx4,jumpTime,Sum90.velSide{jump})
        vsSum(jump,:) = Sum90.velSide{jump}';
        cueDiff_sum(jump,:) = Sum90.cueDiff{jump}';
    end
   
    
    jumpNum = jump;
    for jump = 1:size(Summ90,1)
        %plot(bx1, jumpTime,Sum90.cueAngle{jump})
        %plot(bx2,jumpTime,Sum90.velFor{jump})
        VmSum(jump + jumpNum,:) = Summ90.jump_vm{jump};
        vfSum(jump + jumpNum,:) = Summ90.velFor{jump}'; 
        %plot(bx3,jumpTime,Sum90.velYaw{jump})
        vySum(jump + jumpNum,:) = Summ90.velYaw{jump}';
        %plot(bx4,jumpTime,Sum90.velSide{jump})
        vsSum(jump + jumpNum,:) = Summ90.velSide{jump}';
        cueDiff_sum(jump + jumpNum,:) = Summ90.cueDiff{jump}';
    end
    
    SEM_Vm = std(VmSum,[],1,'omitnan') / sqrt(size(VmSum,1));
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vs = std(vsSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    nvfSum = zeros(size(nSum90,1),size(jump_vf,1));
    nvySum = zeros(size(nSum90,1),size(jump_vf,1));
    nvsSum = zeros(size(nSum90,1),size(jump_vf,1));
    nVmSum = zeros(size(nSum90,1),size(jump_vf,1));
    ncueDiff_sum = zeros(size(nSum90,1),size(cueDiff,1));
    for jump = 1:size(nSum90,1)  
        %plot(bx1, jumpTime,Sum90.cueAngle{jump})
        %plot(bx2,jumpTime,Sum90.velFor{jump})
        nVmSum(jump,:) = nSum90.jump_vm{jump};
        nvfSum(jump,:) = nSum90.velFor{jump}'; 
        %plot(bx3,jumpTime,Sum90.velYaw{jump})
        nvySum(jump,:) = nSum90.velYaw{jump}';
        %plot(bx4,jumpTime,Sum90.velSide{jump})
        nvsSum(jump,:) = nSum90.velSide{jump}';
        ncueDiff_sum(jump,:) = nSum90.cueDiff{jump}';
    end
    
    jumpNum = jump;
    for jump = 1:size(nSumm90,1)
        %plot(bx1, jumpTime,Sum90.cueAngle{jump})
        %plot(bx2,jumpTime,Sum90.velFor{jump})
        nVmSum(jump + jumpNum,:) = nSumm90.jump_vm{jump};
        nvfSum(jump + jumpNum,:) = nSumm90.velFor{jump}'; 
        %plot(bx3,jumpTime,Sum90.velYaw{jump})
        nvySum(jump + jumpNum,:) = nSumm90.velYaw{jump}';
        %plot(bx4,jumpTime,Sum90.velSide{jump})
        nvsSum(jump + jumpNum,:) = nSumm90.velSide{jump}';
        ncueDiff_sum(jump + jumpNum,:) = nSumm90.cueDiff{jump}';
    end
    
    nSEM_Vm = std(nVmSum,[],1,'omitnan') / sqrt(size(nVmSum,1));
    nSEM_vf = std(nvfSum,[],1,'omitnan') / sqrt(size(nvfSum,1));
    nSEM_vy = std(nvySum,[],1,'omitnan') / sqrt(size(nvfSum,1));
    nSEM_vs = std(nvsSum,[],1,'omitnan') / sqrt(size(nvfSum,1));
    
%     figure();
%     sgtitle('90 jumps')
%     set(gcf,'color','w')
%     set(gcf,'renderer','painters')
%     bx1 = subplot(1,1,1); 
%     hold on
    

    keepIndex = ~isnan(nSEM_Vm);
    SEMhigh = [smoothdata(mean(abs(nVmSum(:,keepIndex)),'omitnan'),'gaussian',1) + nSEM_Vm(keepIndex)]; 
    SEMlow = [smoothdata(mean(abs(nVmSum(:,keepIndex)),'omitnan'),'gaussian',1) - nSEM_Vm(keepIndex)];
   %plot(ax1,jumpTime,VmSum,'color',[0.5,0.5,0.5],'LineWidth',0.5)
     plot(bx1,jumpTime,(mean(abs(nVmSum),1)),'r','LineWidth',1.5)
   % plot(bx1,jumpTime,((nVmSum)),'color',[0,0,0,0.2])
    %plot(bx1,jumpTime,var(VmSum,1),'k','LineWidth',1.5)
    patch(bx1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0,0],'FaceAlpha',0.25,'EdgeColor','none')
    xline(bx1,0)
    xlim(bx1,[-plotWin,plotWin])
    ylabel(bx1,'|Z(Vm)|')
    xlabel(bx1,'Time (s)')
%     
%     figure();
%     sgtitle('90 jumps')
%     set(gcf,'color','w')
%     set(gcf,'renderer','painters')
%     bx1 = subplot(1,1,1); 
%     hold on
    
    keepIndex = ~isnan(SEM_Vm);
    SEMhigh = [smoothdata(mean(abs(VmSum(:,keepIndex)),'omitnan'),'gaussian',1) + SEM_Vm(keepIndex)]; 
    SEMlow = [smoothdata(mean(abs(VmSum(:,keepIndex)),'omitnan'),'gaussian',1) - SEM_Vm(keepIndex)];
%     %plot(ax1,jumpTime,VmSum,'color',[0.5,0.5,0.5],'LineWidth',0.5)
     plot(bx1,jumpTime,(mean(abs(VmSum),1)),'k','LineWidth',1.5)
    % plot(bx1,jumpTime,((VmSum)),'color',[0,0,0,0.1],'LineWidth',0.5)
    %plot(bx1,jumpTime,var(VmSum,1),'k','LineWidth',1.5)
    patch(bx1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xline(bx1,0)
    xlim(bx1,[-plotWin,plotWin])
    ylabel(bx1,'Z(Vm)')
    xlabel(bx1,'Time (s)')
        
    %xlim(ax1,[min(jumpTime),max(jumpTime)])
    keepIndex = ~isnan(SEM_vf);
    SEMhigh = [mean(vfSum(:,keepIndex),'omitnan') + SEM_vf(keepIndex)]; 
    SEMlow = [mean(vfSum(:,keepIndex),'omitnan') - SEM_vf(keepIndex)];
    %plot(ax2,jumpTime,vfSum,'color',[0.5,0.5,0.5],'LineWidth',0.5)
    plot(bx3,jumpTime,mean(vfSum,1),'k','LineWidth',1)
    patch(bx3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xline(bx3,0)
    xlim(bx3,[-plotWin,plotWin])
    ylabel(bx3,'forward velocity (mm/s)')
    xlabel(bx3,'Time (s)')
    
    %xlim(ax1,[min(jumpTime),max(jumpTime)])
    keepIndex = ~isnan(nSEM_vf);
    SEMhigh = [mean(nvfSum(:,keepIndex),'omitnan') + nSEM_vf(keepIndex)]; 
    SEMlow = [mean(nvfSum(:,keepIndex),'omitnan') - nSEM_vf(keepIndex)];
    %plot(ax2,jumpTime,vfSum,'color',[0.5,0.5,0.5],'LineWidth',0.5)
    plot(bx3,jumpTime,mean(nvfSum,1),'r','LineWidth',1.5)
    patch(bx3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0,0],'FaceAlpha',0.25,'EdgeColor','none')
    xline(bx3,0)
    xlim(bx3,[-plotWin,plotWin])
    ylabel(bx3,'forward velocity (mm/s)')
    xlabel(bx3,'Time (s)')
    
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [mean(abs(vySum(:,keepIndex)),'omitnan') + SEM_vy(keepIndex)]; 
    SEMlow = [mean(abs(vySum(:,keepIndex)),'omitnan') - SEM_vy(keepIndex)];
    %plot(ax3,jumpTime,abs(vySum),'color',[0.5,0.5,0.5],'LineWidth',0.5)
    plot(bx2,jumpTime,mean(abs(vySum),1),'k','LineWidth',1.5)
    patch(bx2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xline(bx2,0)
    xlim(bx2,[-plotWin,plotWin])
    ylabel(bx2,'rotational speed (deg/s)')
    xlabel(bx2,'Time (s)')
    
    keepIndex = ~isnan(nSEM_vy);
    SEMhigh = [mean(abs(nvySum(:,keepIndex)),'omitnan') + nSEM_vy(keepIndex)]; 
    SEMlow = [mean(abs(nvySum(:,keepIndex)),'omitnan') - nSEM_vy(keepIndex)];
    %plot(ax3,jumpTime,abs(vySum),'color',[0.5,0.5,0.5],'LineWidth',0.5)
    plot(bx2,jumpTime,mean(abs(nvySum),1),'r','LineWidth',1.5)
    patch(bx2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0,0],'FaceAlpha',0.25,'EdgeColor','none')
    xline(bx2,0)
    xlim(bx2,[-plotWin,plotWin])
    ylabel(bx2,'rotational speed (deg/s)')
    xlabel(bx2,'Time (s)')

linkaxes([bx1,bx2,bx3],'x')

figure; 
sgtitle('90 jumps')
set(gcf,'color','w')
set(gcf,'renderer','painters')
yyaxis left
hold on
keepIndex = ~isnan(SEM_Vm);
SEMhigh = [mean(abs(VmSum(:,keepIndex)),'omitnan') + SEM_Vm(keepIndex)]; 
SEMlow = [mean(abs(VmSum(:,keepIndex)),'omitnan') - SEM_Vm(keepIndex)];
patch([jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,1],'FaceAlpha',0.1,'EdgeColor','none')
plot(jumpTime,(mean(abs(VmSum),1)),'b','LineWidth',1.5)
xline(0)
%ylim([min(SEMlow),max(SEMhigh)])
%ylim([0.4,2.4])
ylabel('|Z(Vm)|')
yyaxis right
hold on
keepIndex = ~isnan(SEM_vy);
SEMhigh = [mean(abs(vySum(:,keepIndex)),'omitnan') + SEM_vy(keepIndex)]; 
SEMlow = [mean(abs(vySum(:,keepIndex)),'omitnan') - SEM_vy(keepIndex)];
patch([jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[1,0,0],'FaceAlpha',0.1,'EdgeColor','none')
plot(jumpTime,mean(abs(vySum),1),'r','LineWidth',1.5)
ylim([min(SEMlow),max(SEMhigh)])
ylabel('|Yaw|')
xlim([-0.2,1])




figure;
set(gcf,'renderer','painters','color','w')
subplot(2,1,1);
sgtitle('90 jumps')
plot(jumpTime,VmSum,'k')
xlim([-plotWin,plotWin])
box off
subplot(2,1,2)
plot(jumpTime,nVmSum,'r')
xlim([-plotWin,plotWin])
box off
    


 

%% uncomment to plot -90 sep from 90
    vfSum = zeros(size(Summ90,1),size(jump_vf,1));
%     vySum = zeros(size(Summ90,1),size(jump_vf,1));
%     vsSum = zeros(size(Summ90,1),size(jump_vf,1)); 
%     VmSum = zeros(size(Summ90,1),size(jump_vf,1)); 
%     figure();
%     sgtitle('-90 jumps')
%     set(gcf,'color','w')
%     set(gcf,'renderer','painters')
%     %cx1 = subplot(4,1,1);
%     %hold on
%     cx1 = subplot(3,1,1);
%     cx2 = subplot(3,1,2);
%     hold on
%     cx3 = subplot(3,1,3); 
%     hold on
% %     cx4 = subplot(3,1,3); 
% %     hold on
%     for jump = 1:size(Summ90,1)  
%         %plot(cx1, jumpTime,Summ90.cueAngle{jump})
%         %plot(cx2,jumpTime,Summ90.velFor{jump})
%         VmSum(jump,:) = Summ90.jump_vm{jump};
%         vfSum(jump,:) = Summ90.velFor{jump}'; 
%         %plot(cx3,jumpTime,Summ90.velYaw{jump})
%         vySum(jump,:) = Summ90.velYaw{jump}';
%         %plot(cx4,jumpTime,Summ90.velSide{jump})
%         vsSum(jump,:) = Summ90.velSide{jump}';
%     end
%     SEM_Vm = std(VmSum,[],1,'omitnan') / sqrt(size(VmSum,1));
%     SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
%     SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
%     SEM_vs = std(vsSum,[],1,'omitnan') / sqrt(size(vfSum,1));
%     
%     keepIndex = ~isnan(SEM_Vm);
%     SEMhigh = [smoothdata(mean(abs(VmSum(:,keepIndex)),'omitnan'),'gaussian',1) + SEM_Vm(keepIndex)]; 
%     SEMlow = [smoothdata(mean(abs(VmSum(:,keepIndex)),'omitnan'),'gaussian',1) - SEM_Vm(keepIndex)];
%     %plot(ax1,jumpTime,VmSum,'color',[0.5,0.5,0.5],'LineWidth',0.5)
%     plot(cx1,jumpTime,mean(var(VmSum),1),'k','LineWidth',1.5)
%     patch(cx1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%     xline(cx1,0)
%     xlim(cx1,[min(jumpTime),max(jumpTime)])
%     ylabel(cx1,'Vm')
%     xlabel(cx1,'Time (s)')
%     
%     
%     %xlim(ax1,[min(jumpTime),max(jumpTime)])
%     keepIndex = ~isnan(SEM_vf);
%     SEMhigh = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',150) + SEM_vf(keepIndex)]; 
%     SEMlow = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',150) - SEM_vf(keepIndex)];
%     %plot(ax2,jumpTime,vfSum,'color',[0.5,0.5,0.5],'LineWidth',0.5)
%     plot(cx3,jumpTime,mean(vfSum,1),'k','LineWidth',1.5)
%     patch(cx3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%     xline(cx3,0)
%     xlim(cx3,[min(jumpTime),max(jumpTime)])
%     ylabel(cx3,'forward velocity (mm/s)')
%     xlabel(cx3,'Time (s)')
%     
%     keepIndex = ~isnan(SEM_vy);
%     SEMhigh = [smoothdata(mean((vySum(:,keepIndex)),'omitnan'),'gaussian',150) + SEM_vy(keepIndex)]; 
%     SEMlow = [smoothdata(mean((vySum(:,keepIndex)),'omitnan'),'gaussian',150) - SEM_vy(keepIndex)];
%     %plot(ax3,jumpTime,abs(vySum),'color',[0.5,0.5,0.5],'LineWidth',0.5)
%     plot(cx2,jumpTime,mean((vySum),1),'k','LineWidth',1.5)
%     patch(cx2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%     xline(cx2,0)
%     xlim(cx2,[min(jumpTime),max(jumpTime)])
%     ylabel(cx2,'rotational speed (deg/s)')
%     xlabel(cx2,'Time (s)')
% 
% linkaxes([cx1,cx2,cx3],'x')
%     
% %     keepIndex = ~isnan(SEM_vs);
% %     SEMhigh = [smoothdata(mean(vsSum(:,keepIndex),'omitnan'),'gaussian',75) + SEM_vs(keepIndex)]; 
% %     SEMlow = [smoothdata(mean(vsSum(:,keepIndex),'omitnan'),'gaussian',75) - SEM_vs(keepIndex)];
% %     plot(cx4,jumpTime,smoothdata(mean(vsSum,'omitnan'),'gaussian',75),'k','LineWidth',1.5)
% %     patch(cx4,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
% %     yline(cx4,0,'k')   
% %     ylim(cx4,[-1,1])
% %     xlim(cx4,[min(jumpTime),max(jumpTime)])
% %     linkaxes([cx2,cx3,cx4],'x')

end