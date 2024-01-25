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

        for t = 1:numTrials
             plotHeadingVectors(folder,pro_behaviourData{t}, t, 1000, 1.5,savePlots)
             PFL2_3_basicLinePlots_ephys(folder,pro_behaviourData{t}, pro_trialData{t},savePlots)
            angleBin = 10;  
             Activity_Heatmap_ephys(folder,angleBin,pro_behaviourData{t}, pro_trialData{t},savePlots)
%             pause(0.1)
%             close all
%             basic_xcorr_PFL(folder,t,pro_trialData{t}, pro_behaviourData{t},savePlots)
            bData = pro_behaviourData{t}; 
            tData = pro_trialData{t};
            prefHead = calcPrefHead(bData, tData);
            [vf_coeff, vy_coeff] = PFL2_3_behaviourVSactivity_lineplots_new(prefHead, 0.5, 120, tData, bData) ;
            saveas(gcf,fullfile(folder,'figures','headSep_linePlots.fig'))
            Plot2DTrajectory_ephys(bData,folder,0);

            %pause
        end
        
%     catch 
%         disp([folder,' failed']);
%         failed{count} = folder;
%         count = count + 1; 
%     end
end









 



