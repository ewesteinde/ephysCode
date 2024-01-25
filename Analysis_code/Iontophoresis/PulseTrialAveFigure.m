%% plot ave of trials
folders = get_folders_ephys_behaviour(rootDir, 1); 

folderNum = length(folders);
fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
countFail = 1;
for ff = 1:folderNum
          %% Get folder information
          folder = folders(ff).folder;


        load(fullfile(folder,'trialData.mat'));
        %load(fullfile(folder,'pro_behaviourData.mat'));
        load(fullfile(folder,'trialMeta.mat'));
        
%     if ~any(strcmp('time',fieldnames(processed_trialData{t})))
%         time = [0:1/1000:length(processed_trialData{1, 1}.scaledOutput)/1000 - 1/1000];
%     else
%         time = processed_trialData{t}.time;
%     end

    start = 1;
    last1 = length(trialData);
    trialDataAve1 = zeros(size(trialData{1, 1}.scaledOutput));
    figure; clf;
    h(1) = subplot(4,1,1:3);
    for t = start:last1
        trialDataAve1 = trialDataAve1 + trialData{t}.scaledOutput;
        plot(trialData{t}.time, trialData{t}.scaledOutput)
        hold on
    end
    trialDataAve1 = trialDataAve1/t;
    plot(trialData{t}.time, trialDataAve1, 'k', 'LineWidth',1.5) 
    ylabel('Voltage (mV)')
    xlim([trialData{t}.time(1), trialData{t}.time(end)])

    h(2) = subplot(4,1,4);
    plot(trialData{t}.time, trialData{t}.input,'k')
    ylabel('Ionto Pulse')
    xlabel('Time (s)')
    sgtitle(strcat('Ionto Pulse Trials ',num2str(start),'-',num2str(last1), ' 10mM ATP'))
    xlim([trialData{t}.time(1), trialData{t}.time(end)])
    linkaxes(h,'x')

    %%
    saveas(gcf,fullfile(folder,'figures','IontoPulseAverage.fig'))
end