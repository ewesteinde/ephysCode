function percJumpsCorr_rho(rootDir, PFL2, PFL3)       

     folders = get_folders(rootDir,PFL2,PFL3);
     
     headers = {'folder','trial','jumpSize','jumpIdx','velFor','velSide','velYaw','cueAngle','corrected','timeMov','Z','bumpAmp','region','preJump','postJump'};
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
       try
            if strcmp(folder(end),'.')
                folder = folder(1:end-2); 
            end

             processedData_dir = fullfile(folder,'processed_data');

            % Get data files
            expID = get_expID(folder);
            expList = {expID};

            % Load metadata 
            [expMd, trialMd] = load_metadata(expList, folder);

            % Load imaging data
            roiData = load_roi_data(expList, folder);

            % Load FicTrac data
            [~,~, ~] = load_ft_data(expList, folder, 1, 0);

            % Load panels metadata
            panelsMetadata = load_panels_metadata(expList, folder);

            try
                numTrials = max(size(unique(roiData.trialNum),1),length(trialMd.trialNum)); 
            catch
                numTrials = 1; 
            end
        %%
        for nTrial = 1:numTrials
        %% get processed data
            bump_params = [];
            data_filelist = dir(processedData_dir);
            for files = 1:length(data_filelist)
                if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
                    load(fullfile(processedData_dir,data_filelist(files).name));
                end
            end
            
            load(fullfile(processedData_dir,['zscored_df_f_Trial_00',num2str(nTrial),'.mat']))

            %ftT = ftT_minSmooth; 
            jumpWin = 10;
            
            [~, jump_array, ~] = detect_jumps(ftT, jumpWin, jumpWin,0,0);
            [~, jump_array_rho, ~] = detect_jumps(ftT, 30, 15,0,0);
            
            if PFL3 
                activity = ZData.Z{2} - ZData.Z{1};
            else
                all_roi = [];
                for roi = 1:size(ZData,1)
                    all_roi(roi,:) = ZData.Z{roi};
                end
                activity = mean(all_roi,1,'omitnan')';
            end

                count90 = 1;
                count180 = 1;
                jump90_time = []; 
                jump180_time = []; 
                all_roi = []; 
                
                    jumpIdx = [jump_array(1,1):jump_array(1,3)]';  
                    jumpIdx_rho = [jump_array_rho(1,1):jump_array_rho(1,3)]';
                    
                    if jumpIdx(1) < 1 || jumpIdx_rho(1) < 1
                        jump_array(1,:) = []; 
                        jump_array_rho(1,:) = []; 
                    end

                    for jump = 1:size(jump_array,1)
                        jumpIdx = [jump_array(jump,1):jump_array(jump,3)]'; 
                        if jumpIdx(end) > size(ftT.velFor{1},1)
                            jump_array(jump,:) = []; 
                            jump_array_rho(jump,:) = []; 
                        end
                    end 
                    
                for jump = 1:size(jump_array,1)
                    jumpSize = jump_array(jump,4);
                    jumpIdx = [jump_array(jump,1):jump_array(jump,3)]';  
                    jump_vf = ftT.velFor{1}(jumpIdx);
                    jump_vy = (rad2deg(ftT.velYaw{1}(jumpIdx)));
                    jump_vs = (ftT.velSide{1}(jumpIdx));
                    jump_angle = ftT.cueAngle{1}(jumpIdx); 
                    
                    
                    preJumpIdx = [jump_array(jump,1):jump_array(jump,2)-1];
                    preJumpIdx_rho = [jump_array_rho(jump,1):jump_array_rho(jump,2)-1];
                    
                    timeMov = (sum(abs(ftT.velFor{1}(preJumpIdx_rho)) + abs(ftT.velSide{1}(preJumpIdx_rho)) + abs(deg2rad(ftT.velYaw{1}(preJumpIdx_rho))*4.5) > 1.5))/60; % seconds
                    [rho, ~] = CalculateAverageHeading(ftT,1.5, preJumpIdx_rho);
                    [~, preJumpAngle] = CalculateAverageHeading(ftT,0, preJumpIdx);
                    jump_activity = activity(jumpIdx);
                    
                    %preJumpAngle = ftT.cueAngle{1}(preJumpIdx(end));
                    
                    postJumpIdx = jump_array(jump,2)+1;
                    sizePostJump = size(ftT.cueAngle{1}(postJumpIdx:jump_array(jump,3)));
                    preJump = ones(sizePostJump) * preJumpAngle;
                    cueDiff = angdiff(preJump,deg2rad(ftT.cueAngle{1}(postJumpIdx:jump_array(jump,3)))); 
                    
                     if jumpSize == 180
                        correct = sum(abs(cueDiff) < deg2rad(60),'omitnan'); 
                    else
                        correct = sum(abs(cueDiff) < deg2rad(30),'omitnan');
                    end
                    
                    jump_summary.folder(count) = {folder}; 
                    jump_summary.trial(count) = {nTrial}; 
                    jump_summary.jumpSize(count) = {jumpSize}; 
                    jump_summary.jumpIdx(count) = {jumpIdx};
                    jump_summary.velFor(count) = {jump_vf};
                    jump_summary.velYaw(count) = {jump_vy};
                    jump_summary.velSide(count) = {jump_vs}; 
                    jump_summary.cueAngle(count) = {jump_angle};
                    jump_summary.timeMov(count) = {timeMov}; 
                    jump_summary.rho(count) = {rho};
                    jump_summary.correct(count) = {correct};
                    jump_summary.preJump(count) = {ftT.cueAngle{1}(jump_array(jump,2)-1)}; 
                    jump_summary.postJump(count) = {ftT.cueAngle{1}(jump_array(jump,2)+1)};
                    jump_summary.Z(count) = {jump_activity};

                  
                    if  correct >= 1 && timeMov > 1
                        jump_summary.corrected(count) = {1};
                    elseif correct < 1 && timeMov > 1
                        jump_summary.corrected(count) = {0}; 
                    else
                        jump_summary.corrected(count) = {3};
                    end
                    
                    count = count + 1; 
                    
                end
        end
        catch
            disp([folder, ' failed'])
        end
    end 
    
    SumCorr = jump_summary(cell2mat(jump_summary.corrected) == 1,:); 
    SumnCorr = jump_summary( cell2mat(jump_summary.corrected) == 0,:); 
%%

sampRate = 60; 
lags = [-1:1/sampRate:1];
corrWin = 1; % calc corr over x sec of data before & after jump
jumpI = jumpWin*60 + 1; 
winDiff = jumpWin - corrWin;
corrValuesP = zeros(size(SumCorr,1),length(lags),2); 

for t = 1:length(lags)
    lag = round(lags(t)*sampRate); 
    for jump = 1:size(SumCorr,1)
        % + lag, comparing yaw values w/ activity values delayed by lag
        % - lag, comparing yaw values w/ activity values that preceded by lag
        activity = SumCorr.Z{jump}(jumpI - corrWin*sampRate + lag : jumpI + corrWin*sampRate + lag); 
        if PFL3
            yaw = SumCorr.velYaw{jump}(jumpI - corrWin*sampRate : jumpI + corrWin*sampRate);
        else
            yaw = abs(SumCorr.velYaw{jump}(jumpI - corrWin*sampRate : jumpI + corrWin*sampRate));
        end
    
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
 
%%   
    rhoEdges = [0:0.01:1];
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
    histogram(cell2mat(SumCorr.correct),100); 