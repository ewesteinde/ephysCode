folders = get_folders_ephys(rootDir); 
count = 1; 
for f = 1:size(folders,1)
   % try
        folder = folders(f).folder; 
        if strcmp(folder(end),'.')
            folder = folder(1:end-2); 
        end
        load(fullfile(folder,'trialMeta.mat'));
        
        processedDir = fullfile(folder,'processedData');
        if ~exist(processedDir, 'dir')
            mkdir(processedDir)
        end
        
        if isempty(regexp(trialMeta.fly.flyExp, 'kir'))
            try
                load(fullfile(processedDir,'pro_behaviourData.mat'))
                load(fullfile(processedDir,'pro_trialData.mat'))
            catch
                load(fullfile(folder,'pro_behaviourData.mat'))
                load(fullfile(folder,'pro_trialData.mat'))
                pro_behaviourData = processed_behaviourData; 
                pro_trialData = processed_trialData; 
            end
            numTrials = size(pro_behaviourData,1);

            if numTrials == 1
                writetable(pro_behaviourData{1}, fullfile(processedDir, 'bData_trial_1.csv'))
                writetable(pro_trialData{1}, fullfile(processedDir, 'tData_trial_1.csv'))
            else
                for t = 1:numTrials
                    writetable(pro_behaviourData{t}, fullfile(processedDir, 'bData_trial_',num2str(t),'.csv'))
                    writetable(pro_trialData{t}, fullfile(processedDir, 'tData_trial_',num2str(t),'.csv'))
                end
            end
%         else
%             try
%                 load(fullfile(processedDir,'pro_behaviourData.mat'))
%             catch
%                 load(fullfile(folder,'pro_behaviourData.mat'))
%                 pro_behaviourData = processed_behaviourData; 
%             end
%             numTrials = size(pro_behaviourData,1);
% 
%             if numTrials == 1
%                 writetable(pro_behaviourData{1}, fullfile(folder, 'bData_trial_1.csv'))
%             else
%                 for t = 1:numTrials
%                     writetable(pro_behaviourData{t}, fullfile(folder, 'bData_trial_',num2str(t),'.csv'))
%                 end
%             end
        end
end 