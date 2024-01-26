function PulseTrialAveBehaviourFig(rootDir)
%% find appropriate folders
folders_all = get_folders_ephys_behaviour(rootDir, 1); 
count = 1; 
for f = 1:length(folders_all)
    if regexp(folders_all(f).folder, 'ionto_pCL') 
        folders(count) = folders_all(f);
        count = count + 1;
    end
end


%% Process each folder
folderNum = length(folders);
fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
countFail = 1;
for ff = 1:folderNum
      %% Get folder information
      folder = folders(ff).folder;
      
    load(fullfile(folder,'pro_trialData.mat'));
    load(fullfile(folder,'pro_behaviourData.mat'));

%%
%% plot ave of trials
figure; clf;
start = 1;
last1 = length(processed_trialData);
ScOAve1 = zeros(size(processed_trialData{1, 1}.scaledOutput));
h(6) = subplot(16,1,1:3);
for t = start:last1
    ScOAve1 = ScOAve1 + processed_trialData{t}.scaledOutput;
    plot(processed_trialData{t}.time, processed_trialData{t}.scaledOutput)
    hold on
end
ScOAve1 = ScOAve1/t;
plot(processed_trialData{t}.time, ScOAve1, 'k', 'LineWidth',1.5) 
ylabel('V (mV)')
xlim([processed_trialData{t}.time(1), processed_trialData{t}.time(end)])

%% angle 

start = 1;
last1 = length(processed_behaviourData);
h(1) = subplot(16,1,4:6);
for t = start:last1
    plot(processed_trialData{t}.time, processed_behaviourData{t}.angle)
    hold on
end
ylabel('deg')
xlim([processed_trialData{t}.time(1), processed_trialData{t}.time(end)])

%% for
start = 1;
last1 = length(processed_behaviourData);
velForAve1 = zeros(size(processed_behaviourData{1, 1}.vel_for));
h(2) = subplot(16,1,7:9);
for t = start:last1
    velForAve1 = velForAve1 + processed_behaviourData{t}.vel_for;
    plot(processed_trialData{t}.time, processed_behaviourData{t}.vel_for)
    hold on
end
velForAve1 = velForAve1/t;
plot(processed_trialData{t}.time, velForAve1, 'k', 'LineWidth',1.5) 
ylabel('Vf (mm/sec)')
xlim([processed_trialData{t}.time(1), processed_trialData{t}.time(end)])

%% yaw
start = 1;
last1 = length(processed_behaviourData);
velYawAve1 = zeros(size(processed_behaviourData{1, 1}.vel_yaw));
h(3) = subplot(16,1,10:12);
for t = start:last1
    velYawAve1 = velYawAve1 + abs(processed_behaviourData{t}.vel_yaw);
    plot(processed_trialData{t}.time, abs(processed_behaviourData{t}.vel_yaw))
    hold on
end
velYawAve1 = velYawAve1/t;
plot(processed_trialData{t}.time, velYawAve1, 'k', 'LineWidth',1.5) 
ylabel('abs Vy (deg/sec))')
xlim([processed_trialData{t}.time(1), processed_trialData{t}.time(end)])

%% side

start = 1;
last1 = length(processed_behaviourData);
velSideAve1 = zeros(size(processed_behaviourData{1, 1}.vel_side));
h(4) = subplot(16,1,13:15);
for t = start:last1
    velSideAve1 = velSideAve1 + abs(processed_behaviourData{t}.vel_side);
    plot(processed_trialData{t}.time, abs(processed_behaviourData{t}.vel_side))
    hold on
end
velSideAve1 = velSideAve1/t;
plot(processed_trialData{t}.time, velSideAve1, 'k', 'LineWidth',1.5) 
ylabel('abs Vs (mm/sec)')
xlim([processed_trialData{t}.time(1), processed_trialData{t}.time(end)])

h(5) = subplot(16,1,16);
plot(processed_trialData{t}.time, processed_behaviourData{t}.stim,'k')
ylabel('Ionto Pulse')
xlabel('Time (s)')
sgtitle(strcat('Ionto Pulse Trials ',num2str(start),'-',num2str(last1), ' 10mM ATP'))
xlim([processed_trialData{t}.time(1), processed_trialData{t}.time(end)])
linkaxes(h,'x')

saveas(gcf,fullfile(folder,'figures','expAve.fig'))
end
end

