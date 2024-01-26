folders = get_folders_ephys(rootDir);
savePlots = 1; 
count = 1; 
for f = 1:size(folders,1)
   % try
        folder = folders(f).folder; 
        if strcmp(folder(end),'.')
            folder = folder(1:end-2); 
        end
        
        try
            processedDir = fullfile(folder,'processedData');
            load(fullfile(processedDir,'pro_behaviourData.mat'))
            load(fullfile(processedDir,'pro_trialData.mat'))
        catch
            load(fullfile(folder,'pro_behaviourData.mat'))
            load(fullfile(folder,'pro_trialData.mat'))
            pro_behaviourData = processed_behaviourData; 
            pro_trialData = processed_trialData; 
        end
        numTrials = size(pro_behaviourData,1);

        %% visualize individual trials

         t = 1;
         try
            plotHeadingVectors(folder,pro_behaviourData{t}, t, 1000, 1.5,savePlots)
         catch
             disp('Heading Vector plot failed')
         end
         PFL2_3_basicLinePlots_ephys(folder,pro_behaviourData{t}, pro_trialData{t},savePlots);
         basicLinePlots_n0vel_ephys(folder,pro_behaviourData{t}, pro_trialData{t},savePlots);
         angleBin = 10;  
         Activity_Heatmap_ephys(folder,angleBin,pro_behaviourData{t}, pro_trialData{t},savePlots);
         plot2DTrajectory_activity(folder, savePlots);
         
         [xPos, yPos] = calcCartesianPos(pro_behaviourData{t});

        pause
end
