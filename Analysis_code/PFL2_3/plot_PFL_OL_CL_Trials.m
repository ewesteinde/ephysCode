folders = get_folders_ephys(rootDir);
clean = 1; 
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
    processed_behaviourData = {};
    processed_trialData = {};
    
    %if ~exist(processedDir,'dir')
        load(fullfile(folder,'CL.mat'))
        load(fullfile(folder,'behaviorData.mat'))
        load(fullfile(folder,'CL_idx.mat'));
        CLidx = []; 
        for i = 1:size(CL_startStopIdx,1)
            CLidx = [CLidx, [CL_startStopIdx(i,1):CL_startStopIdx(i,2)]];
        end
        processed_behaviourData{1}.angle = CL.angle';
        processed_behaviourData{1}.vel_for = CL.vel_for';
        processed_behaviourData{1}.vel_yaw = CL.vel_yaw'; 
        processed_behaviourData{1}.vel_side = CL.vel_side';
        processed_behaviourData{1}.time = [1/1000:1/1000:length(CL.vel_for)/1000]';

        processed_trialData{1}.scaledOutput = CL.scaledOutput_down';
        processed_trialData{1}.smooth_Vm = CL.smooth_Vm';
        processed_trialData{1}.fRate_sec = CL.fRate_sec';
        processed_trialData{1}.time = [1/1000:1/1000:length(CL.vel_for)/1000]';
        %processed_trialData{1} = struct2table(processed_trialData{1}); 
        
         fHz = 1000;
         ephysSettings
        [frY,time] = resample_new(behaviorData{1}.frY,fHz,(settings.sampRate));
        frY = frY(round(CLidx));
        [~, jumps,frY_clean] = detect_jumps_ephys(frY, 10, fHz);
        frY = frY_clean;
        processed_behaviourData{1}.frY = frY;
        processed_behaviourData{1}.jumps = jumps;
        %processed_behaviourData{1} = struct2table(processed_behaviourData{1}); 
        mkdir(fullfile(folder,'processedData'))
%     else
%         load(fullfile(processedDir,'pro_behaviourData.mat'))
%         load(fullfile(processedDir,'pro_trialData.mat'))
%     end

        figure_dir = fullfile(folder,'figures');

        if ~exist(figure_dir,'dir')
            mkdir(fullfile(folder,'figures'))
        end

        if iscell(processed_trialData)
        %for trial = 1:length(processed_trialData)
        trial = 1; 

                bData = processed_behaviourData{trial}; 
                tData = processed_trialData{trial};

                WholeTrialFig(bData,tData,figure_dir,trial,savePlot);

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
            %end
        end
        
        if  ~istable(pro_behaviourData{trial})
            pro_behaviourData{trial} = struct2table(pro_behaviourData{trial}); 
        end
        if  ~istable(pro_trialData{trial})
            pro_trialData{trial} = struct2table(pro_trialData{trial}); 
        end

        if savePro
            save(fullfile(processedDir,'pro_behaviourData.mat'),'pro_behaviourData');
            save(fullfile(processedDir,'pro_trialData.mat'),'pro_trialData') ;
        end
    %end
end
