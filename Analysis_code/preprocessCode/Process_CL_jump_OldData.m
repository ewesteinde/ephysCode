function Process_CL_jump_OldData(rootDir)
%% Load in exp data
%close all 
dbstop if error
folders = get_folders_ephys_behaviour(rootDir, 1); 

%% Process each folder
folderNum = length(folders);
fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
countFail = 1;
for ff = 1:folderNum
      %% Get folder information
      folder = folders(ff).folder;
      disp(folder) 
      new = 0; 

    load(fullfile(folder,'trialData.mat'));
    load(fullfile(folder,'trialMeta.mat')); 
    load(fullfile(folder,'behaviorData.mat'));
    load(fullfile(folder,'pro_behaviourData.mat')); 

    if isfield(trialMeta, 'notes')
        disp(trialMeta.notes)
    end
    
    if isfield(trialMeta, 'trialType')
        disp(trialMeta.trialType)
    end
    
    if isfield(trialMeta,'inputR') && isfield(trialMeta,'accessR')
        R = trialMeta.inputR/trialMeta.accessR;
        disp(R)
    end

    if isfield(trialMeta, 'accessREnd')
        Rend = trialMeta.inputREnd/trialMeta.accessREnd;
        disp(Rend)
    end

%% Process Behaviour Data if good cell 
if iscell(trialData) 
    trials = length(trialData);
else
    trials = 1; 
end

for t = 1:trials
ephysSettings

if iscell(behaviorData)
    behaviorData = behaviorData{t};
end
% Interpolate linearly so that original measurements are no longer stepwise, but has smooth transitions

%as of 4/23/21  I only have one smoothing step occuring prior to
%differentiation into velocity. Of all the available methods provided in
%smoothData loess seems to work the best to preserve temporal fidelity
%while reducing noise. Differentiating using the rdiff function rather than
%gradient may produce slightly higher quality traces but requires 5-15
%minutes per 5 minute trial. 

    fHz = 1000;
    [frY,time] = resample_new(behaviorData.frY,fHz,(settings.sampRate));
    [~, jumps,frY_clean] = detect_jumps_ephys(frY, 10, fHz);
    frY = frY_clean;

%     figure(2); clf;
%             k(1) = subplot(4,1,1);
%             plot(time, angle, 'k')
%             ylabel('Pattern Angle')
% 
%             k(2) = subplot(4,1,2);
%             plot(time, vel_for, 'k')
%             ylabel('Vel For mm/s')
% 
%             k(3) = subplot(4,1,3);
%             plot(time, vel_yaw, 'k')
%             ylabel('Vel Yaw deg/s')
% 
%             k(4) = subplot(4,1,4);
%             plot(time, vel_side, 'k')
%             ylabel('Vel Side mm/s')
%             xlabel('Time (s)')
%   
% 
%     figure(3); clf;
%             k(5) = subplot(4,1,1);
%             plot(time, angle, 'k')
%             ylabel('Pattern Angle')
% 
%             k(6) = subplot(4,1,2);
%             plot(time, disp_for, 'k')
%             ylabel('disp For rad/s')
% 
%             k(7) = subplot(4,1,3);
%             plot(time, disp_yaw, 'k')
%             ylabel('disp Yaw rad/s')
% 
%             k(8) = subplot(4,1,4);
%             plot(time, disp_side, 'k')
%             ylabel('disp Side rad/s')
%             xlabel('Time (s)')
%             
%             linkaxes(k,'x');

        if ~istable(processed_behaviourData{t})
            new = 1; 
            processed_behaviourData{t}.frY = frY;
            processed_behaviourData{t}.jumps = jumps;
            processed_behaviourData{t}.vel_for = processed_behaviourData{t}.vel_for';
            processed_behaviourData{t}.vel_yaw = processed_behaviourData{t}.vel_yaw';
            processed_behaviourData{t}.vel_side = processed_behaviourData{t}.vel_side';
            processed_behaviourData{t}.disp_for = processed_behaviourData{t}.disp_for';
            processed_behaviourData{t}.disp_yaw = processed_behaviourData{t}.disp_yaw';
            processed_behaviourData{t}.disp_side = processed_behaviourData{t}.disp_side';
            processed_behaviourData{t}.angle = processed_behaviourData{t}.angle';
            processed_behaviourData{t}.time = processed_behaviourData{t}.time';
            processed_behaviourData{t} = struct2table(processed_behaviourData{t});
        end
        
%     next = input('Next trial? ','s');
%     if ~strcmp(next, 'y')
%         error('fix issue & try again')
%     end
end  

%keep = input('Save? ','s');

if new
    save(fullfile(folder,'pro_behaviourData.mat'),'processed_behaviourData')
end

%% Calculate firing rate
    downsample_Hz = 1000; 
    [processed_trialData] =  calc_fRate(trialData, downsample_Hz, folder);
      
end
end