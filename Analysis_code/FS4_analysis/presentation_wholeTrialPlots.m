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
        
        bData = processed_behaviourData{1};
        tData = processed_trialData{1};
        
%%
        
        xStart = 420;%min(tData.time);
        xEnd = 600;%max(tData.time);         
        
        figure();
        set(gcf,'renderer','painters','color','w','Position', [10 10 1600 600])
        
        a(1) = subplot(4,1,1);
        plot(tData.time, bData.angle,'k')
        xlim([xStart,xEnd])
        ylabel('cue angle')
        
        a(2) = subplot(4,1,2); 
        yyaxis left
        plot(tData.time, tData.scaledOutput,'k')
        ylim([min(tData.scaledOutput),max(tData.scaledOutput)])
        ylabel('Vm (mV)','Color','k')
        set(gca,'YColor','k');
        yyaxis right
        plot(tData.time, bData.vel_for,'r')
        ylabel('vf (mm/s)','Color','r')
        xlabel('Time (s)')
        xlim([xStart,xEnd])
        set(gca,'YColor','r');
        
        a(3) = subplot(4,1,3); 
        yyaxis left
        plot(tData.time, tData.scaledOutput,'k')
        ylim([min(tData.scaledOutput),max(tData.scaledOutput)])
        ylabel('Vm (mV)','Color','k')
        set(gca,'YColor','k');
        yyaxis right
        plot(tData.time, bData.vel_yaw,'r')
        ylabel('vy (deg/s)','Color','r')
        xlabel('Time (s)')
        xlim([xStart,xEnd])
        set(gca,'YColor','r');
        
        a(4) = subplot(4,1,4); 
        yyaxis left
        plot(tData.time, tData.scaledOutput,'k')
        ylim([min(tData.scaledOutput),max(tData.scaledOutput)])
        ylabel('Vm (mV)','Color','k')
        set(gca,'YColor','k');
        yyaxis right
        plot(tData.time, bData.vel_side,'r')
        ylabel('vs (mm/s)','Color','r')
        xlabel('Time (s)')
        xlim([xStart,xEnd])
        set(gca,'YColor','r');
        
        linkaxes(a,'x');
        
        
        
end