function [bestTrials] = findBestSetStimIontoTrials(rootDir, overwrite)

    %folders_all = get_folders_ephys_behaviour(rootDir, 1); 
    
    files = dir(rootDir);
    allTopTrialDir = fullfile(rootDir,'topTrials');
    if ~exist(allTopTrialDir, 'dir')
        mkdir(allTopTrialDir)
    end
    folders = files([files.isdir]);
    dates = folders(3:end); % gets rid of . & .. 
    for d = 1:length(dates) % iterate through every date
        dateDir = fullfile(rootDir, dates(d).name);
        date = dateDir(end-7:end);
        datefiles = dir(dateDir); 
        dateFolders = datefiles([datefiles.isdir]); 
        flies = dateFolders(3:end);  % gets rid of . & .. 
        for f = 1:length(flies) % for every date iterate through every fly
            count = 1; 
            setStimTrials = {}; 
            flyDir = fullfile(dateDir,flies(f).name); 
            fly = flyDir(end-4:end);
            alltrials = get_folders_ephys_behaviour(flyDir, overwrite);
%             figure;
%             ax1 = subplot(2,2
            for t = 1:length(alltrials)
                if ~isempty(regexp(alltrials(t).folder, 'ionto_sCL'))
                    setStimTrials{count} = alltrials(t).folder;
                    count = count + 1;
                end
            end
            
            if ~isempty(setStimTrials)
                for trial = 1:length(setStimTrials)
                    trialDir = setStimTrials{trial}(1:end-2); 
                    load(fullfile(trialDir,'pro_trialData.mat'));
                    load(fullfile(trialDir,'pro_behaviourData.mat'));
                    vf = processed_behaviourData{1}.vel_for; 
                    vs = processed_behaviourData{1}.vel_side;
                    vy = processed_behaviourData{1}.vel_yaw / 360 * (2*pi*4.5); % convert from deg to mm
                    cueAngle = processed_behaviourData{1}.angle;
                    totalMov = vf + abs(vs) + abs(vy); 
                    threshold = 1; 
                    perMov = sum(totalMov > threshold)/length(totalMov) * 100; 
                    ScO = processed_trialData{1}.scaledOutput;
                    
                    figure;
                    h(1) = subplot(3,1,1);
                    plot(cueAngle)
                    title([erase(trialDir,rootDir),' setStimNum: ',num2str(trial)],'Interpreter', 'none')
                    h(2) = subplot(3,1,2); 
                    plot(totalMov)
                    title(['Percentage Movement: ',num2str(perMov)])
                    h(3) = subplot(3,1,3); 
                    histogram(vf)
                    title('vf Distribution')
                    
                     % look at histogram of vf distribution
                     % look at percent moving
                     % look at raw voltage
                end
            end
            
            if isempty(setStimTrials)
                topTrialnum = 0;
            else
               topTrialnum = input('Top trial number '); 
            end
               if topTrialnum ~= 0 
                  topTrial = setStimTrials{topTrialnum}(1:end-2);
                  topTrialFlyDir = fullfile(allTopTrialDir,date,fly);
                    if ~exist(topTrialFlyDir, 'dir')
                        mkdir(topTrialFlyDir)
                    end
                  new_topTrial = strrep(topTrial,'\','/'); % replace backslashes b/c matlab uses them as escape characters
                  fullPathTxt = fopen(fullfile(topTrialFlyDir,'pathName.txt'),'w');
                  fprintf(fullPathTxt,new_topTrial);
                  fclose(fullPathTxt);
                  copyfile(topTrial, topTrialFlyDir)
               end
               close all 
        end
    end
  
end