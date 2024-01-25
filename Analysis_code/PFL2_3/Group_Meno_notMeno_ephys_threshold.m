function [MenoDataFly, nMenoData] = Group_Meno_notMeno_ephys_threshold(rootDir,saveData)
varNames = ["Folder","numTrial","Goal","rho","timeMov","Indices"];
varTypes = ["string","double","double","double","double","cell"];
MenoDataFly = table('Size',[1 6],'VariableNames', varNames,'VariableTypes',varTypes);
nMenoData = table('Size',[1 6],'VariableNames', varNames,'VariableTypes',varTypes);


folders = get_folders_ephys(rootDir);

if ~isempty(regexp(rootDir,'CL_OL'))
    window = 30; 
else
    window = 14; 
end
minVel = 1.5; 
highThres = 0.85; 
lowThres = 0.6; 

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
    mCount = 1;
    nCount = 1;
    MenoDataFly = table('Size',[1 6],'VariableNames', varNames,'VariableTypes',varTypes);
    nMenoDataFly = table('Size',[1 6],'VariableNames', varNames,'VariableTypes',varTypes);
  
    Meno_chunks = [];
    not_Meno_chunks = []; 
    Meno_sum = {};
    not_Meno_sum = {};
    bData = pro_behaviourData{nTrial};
    tData = pro_trialData{nTrial};
    [~, ~, Meno_chunks, not_Meno_chunks,bData_interp] = SegmentMenovsNotMeno_ephys(bData,window, minVel,highThres,lowThres,folder);
    
    if ~isempty(Meno_chunks) || ~isempty(not_Meno_chunks)
        Meno_chunks = Meno_chunks';
        not_Meno_chunks = not_Meno_chunks';

        count = 1; 
        for chunk = 1:length(Meno_chunks)
            interpIdx = Meno_chunks{chunk,1};
            interpTime = bData_interp.time(interpIdx);
            idx = find(bData.time >= interpTime(1) & bData.time <= interpTime(end));
            [rho, theta] = CalculateAverageHeading_ephys(bData,minVel, idx);
            if ~isnan(theta)
                Meno_sum{count,1} = idx'; 
                Meno_sum{count,2} = theta; 
                Meno_sum{count,3} = rho;
                count = count + 1; 
            end
        end
        count = 1; 
        for chunk = 1:length(not_Meno_chunks)
            interpIdx = not_Meno_chunks{chunk,1};
            interpTime = bData_interp.time(interpIdx);
            idx = find(bData.time > interpTime(1) & bData.time < interpTime(end));
            [rho, theta] = CalculateAverageHeading_ephys(bData,minVel, idx);
            if ~isnan(theta)
                not_Meno_sum{count,1} = idx'; 
                not_Meno_sum{count,2} = theta; 
                not_Meno_sum{count,3} = rho;
                count = count + 1; 
            end
        end
    

        % collect chunks with sim goals
            if ~isempty(Meno_sum)
                MenoDataFly.Indices{mCount} = [Meno_sum{:,1}];
                [rho, theta] = CalculateAverageHeading_ephys(bData,minVel, MenoDataFly.Indices{mCount});
                speed = sqrt(bData.vel_side(MenoDataFly.Indices{mCount}).^2 + bData.vel_for(MenoDataFly.Indices{mCount}).^2);
                timeMov = sum(speed > 1.5)/1000;
                MenoDataFly.Goal(mCount) = theta; 
                MenoDataFly.rho(mCount) = rho;
                MenoDataFly.timeMov(mCount) = timeMov; 
                MenoDataFly.Folder(mCount) = folder; 
                MenoDataFly.numTrial(mCount) = nTrial;
            end
            
            if ~isempty(not_Meno_sum)
                nMenoDataFly.Indices{nCount} = [not_Meno_sum{:,1}];
                [rho, theta] = CalculateAverageHeading_ephys(bData,minVel, nMenoDataFly.Indices{nCount});
                speed = sqrt(bData.vel_side(nMenoDataFly.Indices{nCount}).^2 + bData.vel_for(nMenoDataFly.Indices{nCount}).^2);
                timeMov = sum(speed > 1.5)/1000;
                nMenoDataFly.Goal(nCount) = theta; 
                nMenoDataFly.rho(nCount) = rho;
                nMenoDataFly.timeMov(nCount) = timeMov;
                nMenoDataFly.Folder(nCount) = folder;
                nMenoDataFly.numTrial(nCount) = nTrial;
            end
        
        if saveData
            save(fullfile(folder,'processedData',['MenoDataFly_trial_',num2str(nTrial),'.mat']),'MenoDataFly');
            save(fullfile(folder,'processedData',['nMenoDataFly_trial_',num2str(nTrial),'.mat']),'nMenoDataFly');
        end
        
    end
    end
    
end

% if saveData
%         save(fullfile(rootDir,'MenoData.mat'),'MenoData')
%         save(fullfile(rootDir,'nMenoData.mat'),'nMenoData')
% end
    
end