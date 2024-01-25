folders = get_folders_ephys_behaviour(rootDir, 1); 

folderNum = length(folders);
fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
countFail = 1;
for ff = 1:folderNum
    
    folder = folders(ff).folder;
    disp(folder)
    
    load(fullfile(folder,'pro_trialData.mat'));
    load(fullfile(folder,'pro_behaviourData.mat'));
    load(fullfile(folder,'trialMeta.mat'));


t = 1;

%
nActivity = processed_trialData{t}.scaledOutput;
figure();clf;      
        
        h(1) = subplot(5,1,1);
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.angle, 'k') 
        ylabel('angle')
        xlim([processed_behaviourData{t}.time(1), processed_behaviourData{t}.time(end)])
       
        
        h(2) = subplot(5,1,2);
        yyaxis left
        plot(processed_behaviourData{t}.time,nActivity, 'k') 
        ylabel('Vm')
        %ylim([-70 -45])
        hold on 
        yyaxis right
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_for, 'r')
        ylim([-(max(processed_behaviourData{t}.vel_for)) max(processed_behaviourData{t}.vel_for)])
        xlim([processed_behaviourData{t}.time(1), processed_behaviourData{t}.time(end)])
        ylabel('Vf mm/sec')
        
        
        
        h(3) = subplot(5,1,3);
        yyaxis left
        plot(processed_behaviourData{t}.time,nActivity, 'k')
        ylabel('Vm')
       % ylim([-70 -45])
        hold on 
        yyaxis right
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_side, 'r')
        ylim([-(max(processed_behaviourData{t}.vel_side)) max(processed_behaviourData{t}.vel_side)])
        ylabel('Vs mm/sec')
        xlim([processed_behaviourData{t}.time(1), processed_behaviourData{t}.time(end)])
        
        h(4) = subplot(5,1,4);
        yyaxis left
        plot(processed_behaviourData{t}.time,nActivity, 'k')
        ylabel('Vm')
        xlim([processed_behaviourData{t}.time(1), processed_behaviourData{t}.time(end)])

        %ylim([-55 -])
        hold on 
        yyaxis right
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_yaw, 'r')
        ylim([(min(processed_behaviourData{t}.vel_yaw)) max(processed_behaviourData{t}.vel_yaw)])
        ylabel('Vy deg/sec')
        
        h(5) = subplot(5,1,5);
        yyaxis left
        plot(processed_behaviourData{t}.time,nActivity, 'k')
        ylabel('Vm')
        xlim([processed_behaviourData{t}.time(1), processed_behaviourData{t}.time(end)])

        %ylim([-55 -])
        hold on 
        yyaxis right
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.stim, 'r')
        ylabel('Ionto Stim')
        xlabel('Time (s)')
        xlim([processed_behaviourData{t}.time(1), processed_behaviourData{t}.time(end)])
        
        
        linkaxes(h,'x');
        
        %%
        saveas(gcf, fullfile(folder,'figures','wholeTrial.fig'))
end