function Process_CL_jump_Data(rootDir)
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

if isempty(regexp(trialMeta.fly.flyExp, 'kir'))
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
processed_behaviourData = {};
if iscell(trialData) 
    trials = length(trialData);
else
    trials = 1; 
end

for t = 1:length(trials)
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
    %[y, t] = resample_new(x, fs_new, fs_old)
    if length(behaviorData.vel_for) == length(behaviorData.frY)
        [vel_for,~] = resample_new(behaviorData.vel_for,fHz,(settings.sampRate));
        [vel_yaw,~] = resample_new(behaviorData.vel_yaw,fHz,(settings.sampRate));
        [vel_side,~] = resample_new(behaviorData.vel_side,fHz,(settings.sampRate));
    else
        vel_for = behaviorData.vel_for; 
        vel_yaw = behaviorData.vel_yaw;
        vel_side = behaviorData.vel_side;
    end
    
    if length(behaviorData.disp_for) == length(behaviorData.frY)
        [disp_for,~] = resample_new(behaviorData.disp_for,fHz,(settings.sampRate));
        [disp_yaw,~] = resample_new(behaviorData.disp_yaw,fHz,(settings.sampRate));
        [disp_side,~] = resample_new(behaviorData.disp_side,fHz,(settings.sampRate));
    else
        disp_for = behaviorData.disp_for; 
        disp_yaw = behaviorData.disp_yaw;
        disp_side = behaviorData.disp_side;
    end
  
    [angle,~] = resample_new(behaviorData.angle,fHz,(settings.sampRate));
    angle = smoothdata(angle,'movmedian',20);
    angle(angle > 180) = 180;
    angle(angle < -180) = -180; 
%     int_angle = int_angle';

    [frY,time] = resample_new(behaviorData.frY,fHz,(settings.sampRate));
    [~, jumps,frY_clean] = detect_jumps_ephys(frY, 10,10, fHz);
    frY = frY_clean;

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

        processed_behaviourData{t} = table(time, disp_for, disp_side, disp_yaw, angle, vel_for, vel_side, vel_yaw, frY, jumps);
        
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
for t = 1:1%length(behaviorData)

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
            ylim([-180,180])
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
             try
                ylim([(min(processed_behaviourData{t}.vel_yaw)) max(processed_behaviourData{t}.vel_yaw)])
            catch
                ylim([-1 1])
             end
            ylabel('Vy deg/sec')


            linkaxes(h,'x');
        end
        
        fig_dir = fullfile(folder,'figures'); 
        if ~isfolder(fig_dir)
            mkdir(fig_dir)
        end
        
        saveas(gcf, fullfile(fig_dir,['whole_trial',num2str(t),'.fig']))
end

else
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
    processed_behaviourData = {};
    if iscell(trialData) 
        trials = length(trialData);
    else
        trials = 1; 
    end

    fig_dir = fullfile(folder,'figures'); 
    if ~isfolder(fig_dir)
        mkdir(fig_dir)
    end
    behaviorData_all = behaviorData;
    for t = 1:trials
    ephysSettings
    

    if iscell(behaviorData_all)
        behaviorData = behaviorData_all{t};
    end
    % Interpolate linearly so that original measurements are no longer stepwise, but has smooth transitions

    %as of 4/23/21  I only have one smoothing step occuring prior to
    %differentiation into velocity. Of all the available methods provided in
    %smoothData loess seems to work the best to preserve temporal fidelity
    %while reducing noise. Differentiating using the rdiff function rather than
    %gradient may produce slightly higher quality traces but requires 5-15
    %minutes per 5 minute trial. 

    fHz = 1000;
    %[y, t] = resample_new(x, fs_new, fs_old)
    if length(behaviorData.vel_for) == length(behaviorData.frY)
        [vel_for,~] = resample_new(behaviorData.vel_for,fHz,(settings.sampRate));
        [vel_yaw,~] = resample_new(behaviorData.vel_yaw,fHz,(settings.sampRate));
        [vel_side,~] = resample_new(behaviorData.vel_side,fHz,(settings.sampRate));
    else
        vel_for = behaviorData.vel_for; 
        vel_yaw = behaviorData.vel_yaw;
        vel_side = behaviorData.vel_side;
    end
    
    if length(behaviorData.disp_for) == length(behaviorData.frY)
        [disp_for,~] = resample_new(behaviorData.disp_for,fHz,(settings.sampRate));
        [disp_yaw,~] = resample_new(behaviorData.disp_yaw,fHz,(settings.sampRate));
        [disp_side,~] = resample_new(behaviorData.disp_side,fHz,(settings.sampRate));
    else
        disp_for = behaviorData.disp_for; 
        disp_yaw = behaviorData.disp_yaw;
        disp_side = behaviorData.disp_side;
    end
  
    [angle,~] = resample_new(behaviorData.angle,fHz,(settings.sampRate));
    angle = smoothdata(angle,'movmedian',20);
    angle(angle > 180) = 180;
    angle(angle < -180) = -180; 
%     int_angle = int_angle';

    [frY,time] = resample_new(behaviorData.frY,fHz,(settings.sampRate));
    [~, jumps,frY_clean] = detect_jumps_ephys(frY, 10, fHz);
    frY = frY_clean;

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
            
            linkaxes(k,'x');
            
    saveas(gcf, fullfile(fig_dir,['whole_trial',num2str(t),'.fig']))
  

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

        processed_behaviourData{t} = table(time, disp_for, disp_side, disp_yaw, angle, vel_for, vel_side, vel_yaw, frY, jumps);
        
    %     next = input('Next trial? ','s');
    %     if ~strcmp(next, 'y')
    %         error('fix issue & try again')
    %     end
    end  

    %keep = input('Save? ','s');

    %if strcmp(keep, 'y')
        save(fullfile(folder,'pro_behaviourData.mat'),'processed_behaviourData')
    %end
    
end

end
end