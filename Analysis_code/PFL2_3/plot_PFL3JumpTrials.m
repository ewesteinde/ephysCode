folders = get_folders_ephys(rootDir);
clean = 0; 
savePro = 1;
savePlot = 1; 

dbstop if error

for f = 1:length(folders)
    folder = folders(f).folder; 
    pro_behaviourData = {};
    pro_trialData = {};
    
    if strcmp(folder(end),'.')
        folder = folder(1:end-2);
    end
    
    processedDir = fullfile(folder,'processedData');
    
    %if ~exist(processedDir,'dir')
        load(fullfile(folder,'pro_behaviourData.mat'))
        load(fullfile(folder,'pro_trialData.mat'))
        mkdir(fullfile(folder,'processedData'))
%     else
%         load(fullfile(processedDir,'pro_behaviourData.mat'))
%         load(fullfile(processedDir,'pro_trialData.mat'))
%         processed_behaviourData = pro_behaviourData; 
%         processed_trialData = pro_trialData; 
%     end

        figure_dir = fullfile(folder,'figures');

        if ~exist(figure_dir,'dir')
            mkdir(fullfile(folder,'figures'))
        end

        if iscell(processed_trialData)
            for trial = 1:length(processed_trialData)

                bData = processed_behaviourData{trial}; 
                tData = processed_trialData{trial};

                WholeTrialFig(bData,tData,figure_dir,trial,savePlot)

                if clean
                    cutFile = input('cut trial? 0/1 ');
                    if cutFile 
                        cutStart = input('time when to start: ');
                        cutEnd = input('time when to end: ');
                        time = bData.time;

                        [startVal, startIdx] = min(abs(time-cutStart)); 
                        [endVal, endIdx] = min(abs(time-cutEnd)); 

                        fn = fieldnames(bData);
                        for field = 1:numel(fn)
                            bData.(fn{field}) = bData.(fn{field})(startIdx:endIdx);
                        end
                        fn = fieldnames(tData);
                        for field = 1:numel(fn)
                            tData.(fn{field}) = tData.(fn{field})(startIdx:endIdx);
                        end                   
                    end
                end

                Plot2DTrajectory_ephys(folder, bData, trial, 1, savePlot);
                plotHeadingVectors(folder, bData, trial, 1000,1.5,savePlot);

                pro_behaviourData{trial} = bData;
                pro_trialData{trial} = tData;
            end
        end

        if savePro
            save(fullfile(processedDir,'pro_behaviourData.mat'),'pro_behaviourData')
            save(fullfile(processedDir,'pro_trialData.mat'),'pro_trialData') 
        end
    %end
end
