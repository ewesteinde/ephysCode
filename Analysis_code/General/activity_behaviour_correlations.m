folders = get_folders_ephys_behaviour(rootDir, 1); 

%% Process each folder
folderNum = length(folders);
fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
countFail = 1;
maxLag = 2;
sampRate = 1000;
lags = [-maxLag:0.01:maxLag];

%%
minVel = 0;
corrValuesP = zeros(length(folders),length(lags),3);
c=1;
f = figure();
for ff = 1:folderNum
      %try
      %% Get folder information
        folder = folders(ff).folder;
        disp(folder) 
        load(fullfile(folder,'pro_behaviourData.mat'))
        load(fullfile(folder,'pro_trialData.mat'))
        
        bData = processed_behaviourData{1}; 
        tData = processed_trialData{1};
        total_mov_mm = abs(bData.vel_for) + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw))*4.5;
        no0vel_idx = find(total_mov_mm > minVel);
        total_mov_mm = total_mov_mm(no0vel_idx);
        vm =  tData.smoothVm(no0vel_idx);        
        %vm = tData.fRate_sec(no0vel_idx);  
        vf = (bData.vel_for(no0vel_idx));
        vy = (bData.vel_yaw(no0vel_idx));
        vs = (bData.vel_side(no0vel_idx));
        
        vf_cut = vf(maxLag*sampRate:end - maxLag*sampRate);
        vy_cut = vy(maxLag*sampRate:end - maxLag*sampRate);
        vs_cut = vs(maxLag*sampRate:end - maxLag*sampRate);
        total_mov_mm_cut = total_mov_mm(maxLag*sampRate:end - maxLag*sampRate);

        for lagNum = 1:length(lags)
            lag = lags(lagNum);
            vmLag = circshift(vm, -round(lag * sampRate));
            vmLag_cut = vmLag(maxLag*sampRate:end - maxLag*sampRate); 
            

            [R,~] = corr(vmLag_cut, vf_cut,'Type','Pearson');        
            corrValuesP(ff,lagNum,1) = R; 
            [R,~] = corr(vmLag_cut, vy_cut,'Type','Pearson');        
            corrValuesP(ff,lagNum,2) = R; 
            [R,~] = corr(vmLag_cut, vs_cut,'Type','Pearson');        
            corrValuesP(ff,lagNum,3) = R; 
            [R,~] = corr(vmLag_cut, total_mov_mm_cut,'Type','Pearson');        
            corrValuesP(ff,lagNum,4) = R; 

        end 
    subplot(folderNum,4,c)
    plot(lags, corrValuesP(ff,:,1))
    c = c+1;
    subplot(folderNum,4,c)
    plot(lags, corrValuesP(ff,:,2))
    c = c+1;
    subplot(folderNum,4,c)
    plot(lags, corrValuesP(ff,:,3))
    c = c+1;
    subplot(folderNum,4,c)
    plot(lags, corrValuesP(ff,:,4))
    c = c+1;
end

figure()
subplot(1,4,1)
plot(lags, mean(corrValuesP(:,:,1),1))
title('vf')
subplot(1,4,2)
plot(lags, mean(corrValuesP(:,:,2),1))
title('vy')
subplot(1,4,3)
plot(lags, mean(corrValuesP(:,:,3),1))
title('vs')
subplot(1,4,4)
plot(lags, mean(corrValuesP(:,:,4),1))
title('total speed')