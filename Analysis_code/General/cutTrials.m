%rootDir = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\newData';

folders = get_folders_ephys_behaviour(rootDir, 1); 

%% Process each folder
folderNum = length(folders);
fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
countFail = 1;
for ff = 1:folderNum
      %% Get folder information
      folder = folders(ff).folder;
    disp(folder) 
    load(fullfile(folder,'pro_behaviourData.mat'))
    load(fullfile(folder,'pro_trialData.mat'))

    figure(55);clf; 
    yyaxis left
    plot(processed_trialData{1}.time, processed_trialData{1}.scaledOutput)
    xlim([min(processed_trialData{1}.time),max(processed_trialData{1}.time)])
    yyaxis right
    plot(processed_trialData{1}.time, processed_behaviourData{1}.angle)
    xlim([min(processed_trialData{1}.time),max(processed_trialData{1}.time)])
    cut = input('cut trial? ','s' );
    
    if strcmp(cut,'y')
        cStart = input('cut start: ');
        cEnd = input('cut end: '); 
        
        processed_trialData{1} = processed_trialData{1}(processed_trialData{1}.time > cStart & processed_trialData{1}.time < cEnd,:);
        processed_behaviourData{1} = processed_behaviourData{1}(processed_behaviourData{1}.time > cStart & processed_behaviourData{1}.time < cEnd,:);
        
        figure(55);clf; 
        yyaxis left
        plot(processed_trialData{1}.time, processed_trialData{1}.scaledOutput)
        xlim([min(processed_trialData{1}.time),max(processed_trialData{1}.time)])
        yyaxis right
        plot(processed_trialData{1}.time, processed_behaviourData{1}.angle)
        xlim([min(processed_trialData{1}.time),max(processed_trialData{1}.time)])
    end
    
    detrendTrial = input('detrend the trial? ','s');
    if strcmp(detrendTrial,'y')
        try    
            dt_smooth_Vm = detrend(processed_trialData{1}.smooth_Vm);
            baseValue = median(processed_trialData{1}.smooth_Vm(1:(10*100)));

            figure(66);clf;
            subplot(3,1,1)
            plot(processed_trialData{1}.smooth_Vm)
            subplot(3,1,2)
            plot(dt_smooth_Vm)
            subplot(3,1,3)
            plot(dt_smooth_Vm + baseValue)

            processed_trialData{1}.smooth_Vm = dt_smooth_Vm + baseValue;
        catch
            dt_smooth_Vm = detrend(processed_trialData{1}.smoothVm);
            baseValue = median(processed_trialData{1}.smoothVm(1:(10*100)));

            figure(66);clf;
            subplot(3,1,1)
            plot(processed_trialData{1}.smoothVm)
            subplot(3,1,2)
            plot(dt_smooth_Vm)
            subplot(3,1,3)
            plot(dt_smooth_Vm + baseValue)
            

            processed_trialData{1}.smoothVm = dt_smooth_Vm + baseValue;
        end
    end
    
    figure(77);clf; 
    subplot(2,1,1);
    yyaxis left
    plot(processed_trialData{1}.time, processed_trialData{1}.scaledOutput)
    xlim([min(processed_trialData{1}.time),max(processed_trialData{1}.time)])
    yyaxis right
    plot(processed_trialData{1}.time, processed_behaviourData{1}.angle)
    xlim([min(processed_trialData{1}.time),max(processed_trialData{1}.time)])
    subplot(2,1,2)
    try
        plot(processed_trialData{1}.time, processed_trialData{1}.smoothVm)
    catch
        plot(processed_trialData{1}.time, processed_trialData{1}.smooth_Vm)
    end
    xlim([min(processed_trialData{1}.time),max(processed_trialData{1}.time)])
    
    good = input('Save trial? ','s');
    if strcmp(good,'y')
        save(fullfile(folder,'pro_behaviourData.mat'),'processed_behaviourData')
        save(fullfile(folder,'pro_trialData.mat'),'processed_trialData')
    else
        disp(['did not update: ',folder])
    end
    
end