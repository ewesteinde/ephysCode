function Plot_Meno_notMeno_ephys(rootDir,savePlots)

folders = get_folders_ephys(rootDir);
minVel = 1.5; 

for f = 1:size(folders,1)
    
    folder = folders(f).folder; 
    if strcmp(folder(end),'.')
        folder = folder(1:end-2); 
    end
    
    processedDir = fullfile(folder,'processedData');
    load(fullfile(processedDir,'pro_behaviourData.mat'))
    load(fullfile(processedDir,'pro_trialData.mat'))
    
    numTrials = length(pro_behaviourData);
        
    for nTrial = 1:numTrials
        nMenoIdx = [];
        MenoIdx = [];
        
        try
            load(fullfile(processedDir,['MenoDataFly_trial_',num2str(nTrial),'.mat']))
            load(fullfile(processedDir,['nMenoDataFly_trial_',num2str(nTrial),'.mat']))
        
            bData = pro_behaviourData{nTrial};
            tData = pro_trialData{nTrial};

            MenoDataFly = MenoDataFly(~ismembertol(MenoDataFly.rho,1,10^-10),:); % gets rid of trials w/ no heading change --> indicates problem
            MenoDataFly = MenoDataFly(MenoDataFly.timeMov > 5,:);

            nMenoDataFly = nMenoDataFly(~ismembertol(nMenoDataFly.rho,1,10^-10),:); % gets rid of trials w/ no heading change --> indicates problem
            nMenoDataFly = nMenoDataFly(nMenoDataFly.timeMov > 5,:);
            try
                nMeno_xcorr = mean_xcorr_PFL(nMenoDataFly, nTrial, tData, bData, 0);
                saveas(nMeno_xcorr,fullfile(folder,'figures',['nMeno_xcorr_trial_',num2str(nTrial),'.fig'])); 

                for chunk = 1:size(nMenoDataFly,1)
                    nMenoIdx = [nMenoIdx, nMenoDataFly.Indices{chunk}];
                end

                nMenobData = bData(nMenoIdx,:);
                nMenotData = tData(nMenoIdx,:);
                lineFig = PFL2_3_basicLinePlots_ephys(folder, nMenobData, nMenotData,0);
                saveas(lineFig,fullfile(folder,'figures',['nMeno_linePlots_trial_',num2str(nTrial),'.fig'])); 
                angleBin = 10;  
                [vf_vs_heatmap,vy_heatmap] = Activity_Heatmap_ephys(folder, angleBin,nMenobData, nMenotData,0);
                saveas(vf_vs_heatmap,fullfile(folder,'figures',['nMeno_vfvsHeatmap_trial_',num2str(nTrial),'.fig']));
                saveas(vy_heatmap,fullfile(folder,'figures',['nMeno_vyHeatmap_trial_',num2str(nTrial),'.fig']));

            catch 
                disp([folder, 'not meno failed'])
            end
        
%% sep plots for each goal 
            try
                    Meno_xcorr = mean_xcorr_PFL(MenoDataFly, nTrial, tData, bData, 0);
                    saveas(Meno_xcorr,fullfile(folder,'figures',['Meno_xcorr_trial_',num2str(nTrial),'.fig']));
                    for chunk = 1:size(MenoDataFly,1)
                       MenoIdx = [MenoIdx, MenoDataFly.Indices{chunk}];
                    end

                    MenobData = bData(MenoIdx,:);
                    MenotData = tData(MenoIdx,:);
                    lineFig = PFL2_3_basicLinePlots_ephys(folder, MenobData, MenotData,0);
                    saveas(lineFig,fullfile(folder,'figures',['Meno_linePlots_trial_',num2str(nTrial),'.fig'])); 
                    angleBin = 10;  
                    [vf_vs_heatmap,vy_heatmap] = Activity_Heatmap_ephys(folder, angleBin,MenobData, MenotData,0);
                    saveas(vf_vs_heatmap,fullfile(folder,'figures',['Meno_vfvsHeatmap_trial_',num2str(nTrial),'.fig']));
                    saveas(vy_heatmap,fullfile(folder,'figures',['Meno_vyHeatmap_trial_',num2str(nTrial),'.fig']));
                    for chunk = 1:size(MenoDataFly,1)
                        MenoIdx = MenoDataFly.Indices{chunk}; 
                        goal =  rad2deg(MenoDataFly.Goal(chunk));
                        MenobData = bData(MenoIdx,:);
                        MenotData = tData(MenoIdx,:);
                        lineFig = PFL2_3_basicLinePlots_ephys(folder, MenobData, MenotData,0);
                        saveas(lineFig,fullfile(folder,'figures',['Meno_linePlots_trial_',num2str(nTrial),'_Goal_',num2str(goal),'.fig'])); 
                        angleBin = 10;  
                        [vf_vs_heatmap,vy_heatmap] = Activity_Heatmap_ephys(folder, angleBin,MenobData, MenotData,0);
                        saveas(vf_vs_heatmap,fullfile(folder,'figures',['Meno_vfvsHeatmap_trial_',num2str(nTrial),'_Goal_',num2str(goal),'.fig']));
                        saveas(vy_heatmap,fullfile(folder,'figures',['Meno_vyHeatmap_trial_',num2str(nTrial),'_Goal_',num2str(goal),'.fig']));
                        Meno_xcorr = basic_xcorr_PFL(folder,nTrial,MenotData,MenobData,0);
                        saveas(Meno_xcorr,fullfile(folder,'figures',['Meno_xcorr_trial_',num2str(nTrial),'_Goal_',num2str(goal),'.fig']));
                    end

            catch 
                disp([folder, 'meno failed'])
            end
        catch
            disp([folder, 'failed'])
        end
    end
    close all
end