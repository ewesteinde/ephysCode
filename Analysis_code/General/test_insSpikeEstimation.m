folders = get_folders_ephys_behaviour(rootDir, 1); 

%% Process each folder
folderNum = length(folders);
fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
countFail = 1;
for ff = 1:folderNum
      %% Get folder information
    folder = folders(ff).folder;
    close all 

    load(fullfile(folder,'pro_trialData.mat'));

    for t = 1:length(processed_trialData)    

        kernWidth = 20; %ms
        smoothSamples = kernWidth*1; % kernWidth * SampRate (in same units as kernWidth)
        raster = processed_trialData{t}.spikes;
        wingSize = 5.5; % Std
        gSamp = [-(wingSize)*smoothSamples:wingSize*smoothSamples];
        %eSamp = [-smoothSamples:smoothSamples];
        %kernal1 = 1/(smoothSamples *
        %sqrt(2*pi)).*exp(-gSamp.^2/(2*smoothSamples.^2)); gaussian kernel 
        kernel = 1/(smoothSamples * sqrt(2)).*exp(-sqrt(2)*abs(gSamp/smoothSamples)); %exponential kernel

        figure();
        plot(kernel);
        psth = conv(raster,kernel,'same');

        disp(sum(kernel))

        % multiply by sample rate (in sec) to get correct units
        fRate_sec = psth .* 1000; 
        figure(); clf;
        time = processed_trialData{t}.time;

                h(1) = subplot(4,1,1);
                plot(processed_trialData{t}.time, processed_trialData{t}.scaledOutput,'k')
                ylabel('Vm')

                h(2) = subplot(4,1,2);
                plot(time,psth)
                ylabel('PSTH')

                h(3) = subplot(4,1,3);
                plot(time,fRate_sec)
                ylabel('fRate_sec')

                h(4) = subplot(4,1,4);
                plot(time,processed_trialData{t}.smoothVm, 'k')
                ylabel('Vm, no spikes')

                xlabel('Time (s)')

                linkaxes(h,'x');

        processed_trialData{t}.psth = psth;
        processed_trialData{t}.fRate_sec = fRate_sec;
    end
save(fullfile(folder,'pro_trialData.mat'),'processed_trialData')
end
    
