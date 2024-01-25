function Process_CL_Data(rootDir)
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

    load(fullfile(folder,'trialData.mat'));
    load(fullfile(folder,'trialMeta.mat')); 
    load(fullfile(folder,'behaviorData.mat'));
    load(fullfile(folder,'rawData.mat')); 

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

if isstruct(behaviorData)
    behaviorData = {behaviorData(1)}; 
elseif istable(behaviorData)
    behaviorData = {behaviorData}; 
end

if isstruct(trialData)
    trialData = {trialData(1)}; 
elseif istable(trialData)
    if any(strcmp('trialTime',fieldnames(trialData)))
       trialData = renamevars(trialData,'trialTime','time');
    end
    trialData = {trialData}; 
end
behaviorData = cell(1,length(trialData));
for t = 1:length(trialData)
ephysSettings
    if iscell(rawData)
        [disp_for, disp_side, disp_yaw, frX, frY, angle, vel_for, vel_side, vel_yaw] = process_data(rawData{t}, settings.panels.barPattern);
    else
        [disp_for, disp_side, disp_yaw, frX, frY, angle, vel_for, vel_side, vel_yaw] = process_data(rawData, settings.panels.barPattern);
    end
    if length(frX) ~= length(disp_for) %fix for stupid resample error
        frX = resample_new(resample_new(frX,30,2e4),2e4,30);
        frY = resample_new(resample_new(frY,30,2e4),2e4,30);
        angle = resample_new(resample_new(angle,30,2e4),2e4,30);
    end
    
    %try
        behaviorData{t} = table(disp_for, disp_side, disp_yaw, frX, frY, angle, vel_for, vel_side, vel_yaw); 
%     catch
%         disp([folder; "didn't make behaviour table"])
%         load(fullfile(folder,'behaviorData.mat'));
%     end
end
save((fullfile(folder,'behaviorData.mat')), 'behaviorData');

processed_behaviourData = {};
for t = 1:length(trialData)
ephysSettings
% Interpolate linearly so that original measurements are no longer stepwise, but has smooth transitions

%as of 4/23/21  I only have one smoothing step occuring prior to
%differentiation into velocity. Of all the available methods provided in
%smoothData loess seems to work the best to preserve temporal fidelity
%while reducing noise. Differentiating using the rdiff function rather than
%gradient may produce slightly higher quality traces but requires 5-15
%minutes per 5 minute trial. 

    fHz = 1000;
    %[y, t] = resample_new(x, fs_new, fs_old)
    [vel_for,time] = resample_new(behaviorData{t}.vel_for,fHz,(settings.sampRate));
    [vel_yaw,~] = resample_new(behaviorData{t}.vel_yaw,fHz,(settings.sampRate));
    [vel_side,~] = resample_new(behaviorData{t}.vel_side,fHz,(settings.sampRate));
    
    [disp_for,~] = resample_new(behaviorData{t}.disp_for,fHz,(settings.sampRate));
    [disp_yaw,~] = resample_new(behaviorData{t}.disp_yaw,fHz,(settings.sampRate));
    [disp_side,~] = resample_new(behaviorData{t}.disp_side,fHz,(settings.sampRate));
    
    if any(strcmp('input',fieldnames(trialData{t})))
        if length(trialData{t}.input) ~= length(behaviorData{t}.disp_for) %fix for stupid resample error
            fields = fieldnames(trialData{t});
            for f = 1:numel(fields)
                trialData{t}.(fields{f}) = resample_new(resample_new(trialData{t}.(fields{f}),30,2e4),2e4,30);
            end
        end
        [stim,~] = resample_new(trialData{t}.input,fHz,(settings.sampRate));
        stim = smoothdata(stim,'movmedian',20); 
        stim(stim < 0.5 ) = 0 ;
        stim(stim > 0.5 ) = 1;
    elseif any(strcmp('stimPattern',fieldnames(trialData{t})))
        if length(trialData{t}.stimPattern) ~= length(behaviorData{t}.disp_for) %fix for stupid resample error
            trialData{t}.stimPattern = resample_new(resample_new(trialData{t}.stimPattern,30,2e4),2e4,30);
        end
        [stim,~] = resample_new(trialData{t}.stimPattern,fHz,(settings.sampRate));
        stim = smoothdata(stim,'movmedian',20);
        stim(stim < 0.5 ) = 0 ;
        stim(stim > 0.5 ) = 1;
    end
    
    [angle,~] = resample_new(behaviorData{t}.angle,fHz,(settings.sampRate));
    angle = smoothdata(angle,'movmedian',20);
    angle(angle > 180) = 180;
    angle(angle < -180) = -180; 
%     int_angle = int_angle';

    figure(2); clf;
            k(1) = subplot(4,1,1);
            plot(time, angle, 'k')
            ylabel('Pattern Angle')

            k(2) = subplot(4,1,2);
            plot(time, vel_for, 'k')
            ylabel('Vel For mm/s')

            k(3) = subplot(4,1,3);
            plot(time, vel_yaw, 'k')
            ylabel('Vel Yaw deg/s')

            k(4) = subplot(4,1,4);
            plot(time, vel_side, 'k')
            ylabel('Vel Side mm/s')
            xlabel('Time (s)')
  

    figure(3); clf;
            k(5) = subplot(4,1,1);
            plot(time, angle, 'k')
            ylabel('Pattern Angle')

            k(6) = subplot(4,1,2);
            plot(time, disp_for, 'k')
            ylabel('disp For rad/s')

            k(7) = subplot(4,1,3);
            plot(time, disp_yaw, 'k')
            ylabel('disp Yaw rad/s')

            k(8) = subplot(4,1,4);
            plot(time, disp_side, 'k')
            ylabel('disp Side rad/s')
            xlabel('Time (s)')
            
            linkaxes(k,'x');


    if exist('stim','var')
        processed_behaviourData{t} = table(time, disp_for, disp_side, disp_yaw, angle, vel_for, vel_side, vel_yaw, stim);
        
    else
        processed_behaviourData{t} = table(time, disp_for, disp_side, disp_yaw, angle, vel_for, vel_side, vel_yaw); 
    end
    
%     next = input('Next trial? ','s');
%     if ~strcmp(next, 'y')
%         error('fix issue & try again')
%     end
end  

%keep = input('Save? ','s');

%if strcmp(keep, 'y')
    save(fullfile(folder,'pro_behaviourData.mat'),'processed_behaviourData')
%end

%% Calculate firing rate
try
    load(fullfile(folder,'pro_trialData.mat'));
catch 
    downsample_Hz = 1000; 
    [processed_trialData] =  calc_fRate(trialData, downsample_Hz, folder);
end
%% Overlay raw Vm & Vf/Vs 
for t = 1:length(behaviorData)

        nActivity = processed_trialData{t}.scaledOutput;
        
        if max(strcmp('stim',processed_behaviourData{t}.Properties.VariableNames))

            figure(10);clf;      

            h(1) = subplot(5,1,1);
            plot(processed_behaviourData{t}.time, processed_behaviourData{t}.angle, 'k') 
            ylabel('angle')


            h(2) = subplot(5,1,2);
            yyaxis left
            plot(processed_behaviourData{t}.time,nActivity, 'k') 
            ylabel('Vm')
            hold on 
            yyaxis right
            plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_for, 'r')
            ylim([-(max(processed_behaviourData{t}.vel_for)) max(processed_behaviourData{t}.vel_for)])
            ylabel('Vf mm/sec')



            h(3) = subplot(5,1,3);
            yyaxis left
            plot(processed_behaviourData{t}.time,nActivity, 'k')
            ylabel('Vm')
            hold on 
            yyaxis right
            plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_side, 'r')
            ylim([min(processed_behaviourData{t}.vel_side) max(processed_behaviourData{t}.vel_side)])
            ylabel('Vs mm/sec')

            h(4) = subplot(5,1,4);
            yyaxis left
            plot(processed_behaviourData{t}.time,nActivity, 'k')
            ylabel('Vm') 
            %ylim([-65 -45])
            hold on 
            yyaxis right
            plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_yaw, 'r')
            try
                ylim([(min(processed_behaviourData{t}.vel_yaw)) max(processed_behaviourData{t}.vel_yaw)])
            catch
                ylim([-1 1])
            end
            ylabel('Vy deg/sec')
            
            h(5) = subplot(5,1,5);
            yyaxis left
            plot(processed_behaviourData{t}.time,nActivity, 'k')
            ylabel('Vm') 
            %ylim([-65 -45])
            hold on 
            yyaxis right
            plot(processed_behaviourData{t}.time, processed_behaviourData{t}.stim)
            ylabel('stim')


            linkaxes(h,'x');
        else
            figure();clf;      

            h(1) = subplot(4,1,1);
            plot(processed_behaviourData{t}.time, processed_behaviourData{t}.angle, 'k') 
            ylabel('angle')


            h(2) = subplot(4,1,2);
            yyaxis left
            plot(processed_behaviourData{t}.time,nActivity, 'k') 
            ylabel('Vm')
            hold on 
            yyaxis right
            plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_for, 'r')
            ylim([min(processed_behaviourData{t}.vel_for) max(processed_behaviourData{t}.vel_for)])
            ylabel('Vf mm/sec')



            h(3) = subplot(4,1,3);
            yyaxis left
            plot(processed_behaviourData{t}.time,nActivity, 'k')
            ylabel('Vm')
            hold on 
            yyaxis right
            plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_side, 'r')
            ylim([min(processed_behaviourData{t}.vel_side) max(processed_behaviourData{t}.vel_side)])
            ylabel('Vs mm/sec')

            h(4) = subplot(4,1,4);
            yyaxis left
            plot(processed_behaviourData{t}.time,nActivity, 'k')
            ylabel('Vm') 
            %ylim([-65 -45])
            hold on 
            yyaxis right
            plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_yaw, 'r')
            ylim([(min(processed_behaviourData{t}.vel_yaw)) max(processed_behaviourData{t}.vel_yaw)])
            ylabel('Vy deg/sec')


            linkaxes(h,'x');
        end
        
        fig_dir = fullfile(folder,'figures'); 
        if ~isfolder(fig_dir)
            mkdir(fig_dir)
        end
        
        saveas(gcf, fullfile(fig_dir,['whole_trial',num2str(t),'.fig']))
end
        
end
end