function percJumpsCorr_rho_ephys(rootDir)       

     folders1 = [];%get_folders_ephys(rootDir);
     folders2 = get_folders_ephys('Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL3\new_flies');
     
     folders = [folders1; folders2];
     
     headers = {'folder','trial','jumpSize','jumpIdx','velFor','velSide','velYaw','cueAngle','corrected','timeMov','preJump','postJump','jump_vm','cueDiff'};
     jump_summary = cell2table(cell(1,14),'VariableNames',headers); 
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
        jumpWin = 4;

        for t = 1:numTrials

        %[rho_all, ~] = CalcTheta_Rho_ephys(bData,3,30);
        [jump_array_rho, ~, ~] = detect_jumps_ephys(bData.frY, 30,15, 1000);
        [jump_array, ~, ~] = detect_jumps_ephys(bData.frY, jumpWin,jumpWin, 1000);

        total_mov_mm = abs(bData.vel_for + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw)*4.5));
        
        try
            %activity = (tData.smooth_Vm - median(tData.smooth_Vm))/mad(tData.smooth_Vm);
            activity = tData.smooth_Vm;
        catch
            %activity = (tData.smoothVm - median(tData.smoothVm))/mad(tData.smoothVm);
            activity = tData.smoothVm;
        end
        
        %activity = tData.fRate_sec;

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
            jumpIdx_rho = [jump_array_rho(jump,1):jump_array_rho(jump,3)]';
            if jumpIdx(end) > size(bData,1) || jumpIdx_rho(end) > size(bData,1)
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
%                     idx_start = preJumpIdx_rho(1); 
%                     idx_end = preJumpIdx_rho(end); 
                    
%                     padIdx = idx_start - 30*1000:idx_end + 30*1000; 
%                     if padIdx(1) < 1
%                         padIdx = padIdx(padIdx >= 1); 
%                         [rho_all, test_angle] = CalcTheta_Rho_ephys(bData,3,30,padIdx);
%                     elseif padIdx(end) > size(tData,1)
%                         padIdx = padIdx(padIdx <= size(tData,1));
%                         [rho_all, test_angle] = CalcTheta_Rho_ephys(bData,3,30,padIdx);
%                     else
%                         [rho_all, test_angle] = CalcTheta_Rho_ephys(bData,3,30,padIdx);       
%                     end
%                     rho_all = rho_all(padIdx >= idx_start & padIdx <= idx_end)';
%                     test_angle = test_angle(padIdx >= idx_start & padIdx <= idx_end)';
%                     test_angle = sum(test_angle,'omitnan')/length(test_angle); 
                    
                    [rho, ~] = CalculateAverageHeading_ephys(bData,1.5, preJumpIdx_rho);

                    %rho = mean(rho_all,'omitnan'); 
                    
%                     figure(12);clf;
%                     plot(bData.angle(preJumpIdx_rho)); 
%                     ylim([-180 180])
%                     title(['rho: ',num2str(rho)])
                    
                    
                    movIdx = find(abs(bData.vel_for(preJumpIdx)) + abs(bData.vel_side(preJumpIdx)) + abs(deg2rad(bData.vel_yaw(preJumpIdx))*4.5) > 3);
                    movIdx = preJumpIdx(movIdx);
                    timeMov = (sum(abs(bData.vel_for(preJumpIdx_rho)) + abs(bData.vel_side(preJumpIdx_rho)) + abs(deg2rad(bData.vel_yaw(preJumpIdx_rho))*4.5) > 3))/1000; % seconds
                    preJump_base = mean(activity(preJumpIdx));
                    jump_vm = activity(jumpIdx) - preJump_base;
                    %jump_vm = activity(jumpIdx);

                    [~, preJumpAngle] = CalculateAverageHeading_ephys(bData,0, preJumpIdx);
                    
                    
                    %preJumpAngle = ftT.cueAngle{1}(preJumpIdx(end));
                    
                    postJumpIdx = jump_array(jump,2)+ 100;
                    sizePostJump = size(bData.angle(postJumpIdx:jump_array(jump,3)));
                    preJump = ones(sizePostJump) * preJumpAngle;
            
                    
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
                        correct = sum(abs(cueDiff) < deg2rad(75),'omitnan'); 
                    else
                        correct = sum(abs(cueDiff) < deg2rad(40),'omitnan');
                    end 
                    
                    if isempty(correct)
                        pause
                    end
                        
                    
                    jump_summary.folder(count) = {folder}; 
                    jump_summary.trial(count) = {t}; 
                    jump_summary.jumpSize(count) = {jumpSize}; 
                    jump_summary.rho(count) = {rho}; 
                    jump_summary.jumpIdx(count) = {jumpIdx};
                    jump_summary.velFor(count) = {jump_vf};
                    jump_summary.velYaw(count) = {jump_vy};
                    jump_summary.velSide(count) = {jump_vs}; 
                    jump_summary.cueAngle(count) = {jump_angle};
                    jump_summary.timeMov(count) = {timeMov}; 
                    jump_summary.preJump(count) = {bData.angle(jump_array(jump,2)-1)}; 
                    jump_summary.postJump(count) = {bData.angle(jump_array(jump,2)+1)}; 
                    jump_summary.jump_vm(count) = {jump_vm};
                    jump_summary.cueDiff(count) = {cueDiff}; 
                    jump_summary.correct(count) = {correct};

                  
                    if  correct >= 10 && timeMov > 1 && rho > 0.88
                        jump_summary.corrected(count) = {1};
                    elseif correct < 10 && timeMov > 1
                        jump_summary.corrected(count) = {0}; 
                    else
                        jump_summary.corrected(count) = {3};
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
    %%
    SumCorr = jump_summary(cell2mat(jump_summary.corrected) == 1,:);% & abs(cell2mat(jump_summary.jumpSize)) == 90,:); 
    SumnCorr = jump_summary( cell2mat(jump_summary.corrected) == 0,:); 
%     Sum90 = jump_summary(cell2mat(jump_summary.jumpSize) == 90 & cell2mat(jump_summary.corrected) == 1,:);
%     nSum90 = jump_summary(cell2mat(jump_summary.jumpSize) == 90 & cell2mat(jump_summary.corrected) == 0,:);
%     Summ90 = jump_summary(cell2mat(jump_summary.jumpSize) == -90 & cell2mat(jump_summary.corrected) == 1,:);
%     nSumm90 = jump_summary(cell2mat(jump_summary.jumpSize) == -90 & cell2mat(jump_summary.corrected) == 0,:);



%% iterate through jumps @ diff lags & calculate correlation 
sampRate = 1000;
lags = [-1:0.01:1];
corrWin = 1;
corrValuesP = zeros(size(SumCorr,1),length(lags),2); 
jumpI = jumpWin*sampRate + 1;

for t = 1:length(lags)
    lag = round(lags(t)*1000); 
    for jump = 1:size(SumCorr,1)
        % + lag, comparing yaw values w/ activity values delayed by lag
        % - lag, comparing yaw values w/ activity values that preceded by lag
        activity = abs(SumCorr.jump_vm{jump}(jumpI + lag : jumpI + corrWin*sampRate + lag));          
        yaw = abs(SumCorr.velYaw{jump}(jumpI : jumpI + corrWin*sampRate));

%         [R,p] = corr(activity, yaw,'Type','Kendall'); 
%         
%         corrValuesK(jump,t,1) = R; 
%         corrValuesK(jump,t,2) = p; 
        
        [R,p] = corr(activity, yaw,'Type','Pearson'); 
        
        corrValuesP(jump,t,1) = R; 
        corrValuesP(jump,t,2) = p;
        
    end
end

figure();
meanCorr = mean(corrValuesP(:,:,1),1,'omitnan');

SEM_corr = std(corrValuesP(:,:,1),[],1,'omitnan') / sqrt(size(corrValuesP(:,:,1),1));
SEMhigh = [meanCorr + SEM_corr]; 
SEMlow = [meanCorr - SEM_corr];
set(gcf,'renderer','painters','color','w')
%plot(lags,corrValues(:,:,1),'color',[0.75,0.75,0.75],'LineWidth',0.5)
patch([lags fliplr(lags)],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
hold on
plot(lags,meanCorr,'k','LineWidth',1.5)
xlabel('Lag (s)')
ylabel('R')
title('Pearson')
box off


% figure();
% meanCorr = mean(corrValuesP(:,:,2),1,'omitnan');
% 
% SEM_corr = std(corrValuesP(:,:,2),[],1,'omitnan') / sqrt(size(corrValuesP(:,:,2),1));
% SEMhigh = [meanCorr + SEM_corr]; 
% SEMlow = [meanCorr - SEM_corr];
% set(gcf,'renderer','painters','color','w')
% %plot(lags,corrValues(:,:,1),'color',[0.75,0.75,0.75],'LineWidth',0.5)
% patch([lags fliplr(lags)],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
% hold on
% plot(lags,meanCorr,'k','LineWidth',1.5)
% xlabel('Lag (s)')
% ylabel('p-value')
% box off
%%

 %% uncomment to check distribution of corrected jumps  
    rhoEdges = [0:0.02:1];
    figure();
%     histogram(cell2mat(SumnCorr.rho),rhoEdges)
%     hold on
    histogram(cell2mat(SumCorr.rho),rhoEdges)
    xlabel('rho')
    ylabel('num jumps')
    title('Corrected Jumps')
    set(gcf,'color','w','renderer','painters')
    box off
    
    figure();
    histogram(cell2mat(SumnCorr.timeMov),50)
    hold on
    histogram(cell2mat(SumCorr.timeMov),50)

    allrho = [cell2mat(SumnCorr.rho);cell2mat(SumCorr.rho)];
    figure();
    histogram(allrho,rhoEdges)
    
    figure(); 
    histogram(cell2mat(SumCorr.correct),1000); 
    
    
    