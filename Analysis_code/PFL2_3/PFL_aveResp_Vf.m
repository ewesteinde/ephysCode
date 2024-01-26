clear
dbstop if error 
%% jump trials 
myFolder = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL3\CL_jump';
folders_all = dir(myFolder);
folders_all = folders_all(3:end);
%folders = folders_all;


count = 0;
for idx = 1:length(folders_all)
    if contains(folders_all(idx).name,'R')
        count = count + 1; 
        folders(count) = folders_all(idx);
    end
end


lag = -100; 
cellType = 'PFL3';

for exp = 1:length(folders)
	folderName = folders(exp).name;
    [processed_trialData, processed_behaviourData] = loadPFL2_3_CLJump(cellType, folderName);
    
    
 % negative lag means neural activity precedes behaviour
    processed_trialData.fRate_sec = processed_trialData.fRate_sec(1:end+lag);
    processed_trialData.smooth_Vm = processed_trialData.smooth_Vm(1:end+lag);
    processed_behaviourData.angle = processed_behaviourData.angle(1:end+lag);
    
    processed_behaviourData.vel_for = processed_behaviourData.vel_for(abs(lag)+1:end);
    processed_behaviourData.vel_yaw = processed_behaviourData.vel_yaw(abs(lag)+1:end);
    processed_behaviourData.vel_side = processed_behaviourData.vel_side(abs(lag)+1:end);

    minVy = -200;
    maxVy = 200;
    step_Vy = 20;
    
    minVf = -2;
    maxVf = 12;
    step_Vf = 0.5; 

    figure(2);clf;      
        
        h(1) = subplot(2,1,1);
        yyaxis left 
        plot(processed_trialData.scaledOutput_down, 'k') 
        yyaxis right
        plot(processed_behaviourData.angle, 'r') 
        ylabel('angle')
       
        
        h(2) = subplot(2,1,2);
        yyaxis left
        plot(processed_trialData.scaledOutput_down, 'k') 
        ylabel('Vm')
        %ylim([-70 -45])
        hold on 
        yyaxis right
        plot(processed_behaviourData.vel_for, 'r')
        ylim([-(max(processed_behaviourData.vel_for)) max(processed_behaviourData.vel_for)])
        ylabel('Vf mm/sec')

        linkaxes(h,'x');
    
    chop = input('cut out usable portion? y/n ','s');
    if strcmp(chop, 'y')
        iStart = input('starting index '); 
        iEnd = input('ending index '); 
        
        f = fieldnames(processed_behaviourData);
        for k=1:numel(f)
            processed_behaviourData.(f{k}) = processed_behaviourData.(f{k})(iStart:iEnd); 
        end


        f = fieldnames(processed_trialData);
        for k=1:numel(f)
            processed_trialData.(f{k}) = processed_trialData.(f{k})(iStart:iEnd); 
        end
    end
    
%% calculate average VM & FR across Vf values
speedThres = 1.5;
contThres = 1; 
 
[no0Vel_frag, no0Vel_idx] = remove0velocity(speedThres, contThres,  processed_behaviourData);
vy = processed_behaviourData.vel_yaw(no0Vel_idx);
vf = processed_behaviourData.vel_for(no0Vel_idx);

    edges = [minVy:step_Vy:maxVy]; 
    [N, edges, bin] = histcounts(vy, edges);
    tempFR = accumarray((bin+1)', processed_trialData.fRate_sec(no0Vel_idx)', [length(edges) 1]);
    tempVm = accumarray((bin+1)', processed_trialData.smooth_Vm(no0Vel_idx)', [length(edges) 1]);
    mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
    mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
    centers = edges(1:end-1)+diff(edges)/2;

    saveBinsVm_nohead_jump_vy{exp} = [centers' mean_binVm];
    saveBinsFR_nohead_jump_vy{exp} = [centers' mean_binFR];
    
    edges = [minVf:step_Vf:maxVf]; 
    [N, edges, bin] = histcounts(vf, edges);
    tempFR = accumarray((bin+1)', processed_trialData.fRate_sec(no0Vel_idx)', [length(edges) 1]);
    tempVm = accumarray((bin+1)', processed_trialData.smooth_Vm(no0Vel_idx)', [length(edges) 1]);
    mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
    mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
    centers = edges(1:end-1)+diff(edges)/2;

    saveBinsVm_nohead_jump_vf{exp} = [centers' mean_binVm];
    saveBinsFR_nohead_jump_vf{exp} = [centers' mean_binFR];
    
end

%% CL_OL trials

myFolder = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\identified_PFL3\CL_OL';

folders_all = dir(myFolder);
folders_all = folders_all(3:end);
%folders = folders_all;

clear folders
count = 0;
for idx = 1:length(folders_all)
    if contains(folders_all(idx).name,'R')
        count = count + 1; 
        folders(count) = folders_all(idx);
    end
end

cellType = 'PFL3';


for exp = 1:length(folders)
	folderName = folders(exp).name;
    [processed_trialData, processed_behaviourData, fileName] = loadPFL2_3_CLOL(cellType, folderName);
    
    processed_trialData.fRate_sec = processed_trialData.fRate_sec(1:end+lag);
    processed_trialData.smooth_Vm = processed_trialData.smooth_Vm(1:end+lag);
    processed_behaviourData.angle = processed_behaviourData.angle(1:end+lag);
    processed_behaviourData.time = processed_behaviourData.time(1:end+lag);
    processed_trialData.scaledOutput_down = processed_trialData.scaledOutput_down(1:end+lag);
    
    processed_behaviourData.vel_for = processed_behaviourData.vel_for(abs(lag)+1:end);
    processed_behaviourData.vel_yaw = processed_behaviourData.vel_yaw(abs(lag)+1:end);
    processed_behaviourData.vel_side = processed_behaviourData.vel_side(abs(lag)+1:end);
    
    
    minVy = -200;
    maxVy = 200;
    step_Vy = 20;
    
    minVf = -2;
    maxVf = 12;
    step_Vf = 0.5; 

    figure(2);clf;      
        
        h(1) = subplot(2,1,1);
        yyaxis left 
        plot(processed_behaviourData.time,processed_trialData.scaledOutput_down, 'k') 
        yyaxis right
        plot(processed_behaviourData.time,processed_behaviourData.angle, 'r') 
        ylabel('angle')
       
        
        h(2) = subplot(2,1,2);
        yyaxis left
        plot(processed_behaviourData.time,processed_trialData.scaledOutput_down, 'k') 
        ylabel('Vm')
        %ylim([-70 -45])
        hold on 
        yyaxis right
        plot(processed_behaviourData.time,processed_behaviourData.vel_for, 'r')
        ylim([-(max(processed_behaviourData.vel_for)) max(processed_behaviourData.vel_for)])
        ylabel('Vf mm/sec')
        
        linkaxes(h,'x');
        
        figure(3);clf;
        
        g(1) = subplot(2,1,1);
        yyaxis left 
        plot(processed_trialData.scaledOutput_down, 'k') 
        yyaxis right
        plot(processed_behaviourData.angle, 'r') 
        ylabel('angle')
       
        
        g(2) = subplot(2,1,2);
        yyaxis left
        plot(processed_trialData.scaledOutput_down, 'k') 
        ylabel('Vm')
        %ylim([-70 -45])
        hold on 
        yyaxis right
        plot(processed_behaviourData.vel_for, 'r')
        ylim([-(max(processed_behaviourData.vel_for)) max(processed_behaviourData.vel_for)])
        ylabel('Vf mm/sec')

        linkaxes(g,'x');
    
    chop = input('cut out usable portion? y/n ','s');
    if strcmp(chop, 'y')
        iStart = input('starting index '); 
        iEnd = input('ending index '); 
        
        f = fieldnames(processed_behaviourData);
        for k=1:numel(f)
            processed_behaviourData.(f{k}) = processed_behaviourData.(f{k})(iStart:iEnd); 
        end


        f = fieldnames(processed_trialData);
        for k=1:numel(f)
            processed_trialData.(f{k}) = processed_trialData.(f{k})(iStart:iEnd); 
        end
    end
    
%     startOL = input('time 1st OL starts ')*1000; % sec
%     endOL = input('time 1st OL ends ')*1000;
%     lengthCL = 179.99*1000;
%     
%     % save CL & OL segments so I don't have to do this again 
%     [CL, OL, CL_startStopIdx, OL_startStopIdx] = separateCL_OL(startOL, endOL, lengthCL, processed_behaviourData, processed_trialData);
%     fileName_CL = fullfile(fileName, 'CL.mat');
%     fileName_CLidx = fullfile(fileName, 'CL_idx.mat');
%     fileName_OL = fullfile(fileName, 'OL.mat');
%     fileName_OLidx = fullfile(fileName, 'OL_idx.mat');
%     
%     save(fileName_CL,'CL')
%     save(fileName_CLidx,'CL_startStopIdx')
%     save(fileName_OL,'OL')
%     save(fileName_OLidx,'OL_startStopIdx')
    
%% calculate average VM & FR across Vf values
[no0Vel_frag, no0Vel_idx] = remove0velocity(speedThres, contThres,  processed_behaviourData);
vy = processed_behaviourData.vel_yaw(no0Vel_idx);
vf = processed_behaviourData.vel_for(no0Vel_idx);(no0Vel_idx);


    edges = [minVy:step_Vy:maxVy]; 
    [N, edges, bin] = histcounts(vy, edges);
    tempFR = accumarray((bin+1)', processed_trialData.fRate_sec(no0Vel_idx)', [length(edges) 1]);
    tempVm = accumarray((bin+1)', processed_trialData.smooth_Vm(no0Vel_idx)', [length(edges) 1]);
    mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
    mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
    centers = edges(1:end-1)+diff(edges)/2;

    saveBinsVm_nohead_CLOL_vy{exp} = [centers' mean_binVm];
    saveBinsFR_nohead_CLOL_vy{exp} = [centers' mean_binFR];
    
    edges = [minVf:step_Vf:maxVf]; 
    [N, edges, bin] = histcounts(vf, edges);
    tempFR = accumarray((bin+1)', processed_trialData.fRate_sec(no0Vel_idx)', [length(edges) 1]);
    tempVm = accumarray((bin+1)', processed_trialData.smooth_Vm(no0Vel_idx)', [length(edges) 1]);
    mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
    mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
    centers = edges(1:end-1)+diff(edges)/2;

    saveBinsVm_nohead_CLOL_vf{exp} = [centers' mean_binVm];
    saveBinsFR_nohead_CLOL_vf{exp} = [centers' mean_binFR];
    
end

%% Add into array & plot
totalsaveBinsVmVf = [saveBinsVm_nohead_jump_vf saveBinsVm_nohead_CLOL_vf];
totalsaveBinsFRVf = [saveBinsFR_nohead_jump_vf saveBinsFR_nohead_CLOL_vf];

totalsaveBinsVmVy = [saveBinsVm_nohead_jump_vy saveBinsVm_nohead_CLOL_vy];
totalsaveBinsFRVy = [saveBinsFR_nohead_jump_vy saveBinsFR_nohead_CLOL_vy];

edge2plot_vy = [5:16];
edge2plot_vf = [3:15];

sumResp = zeros(length(totalsaveBinsVmVf{1}(edge2plot_vf,1)),1);
count = 0;
figure(1);clf;
for cells = 1:length(totalsaveBinsFRVf)
    if sum(isnan(totalsaveBinsFRVf{cells}(edge2plot_vf,2))) > 0
    else
        count = count + 1;
        plot(totalsaveBinsFRVf{1}(edge2plot_vf,1),normalize(totalsaveBinsFRVf{cells}(edge2plot_vf,2)),'color',[0.75 0.75 0.75],'LineWidth',1)
        sumResp = sumResp + normalize(totalsaveBinsFRVf{cells}(edge2plot_vf,2));
        hold on
    end
end
plot(totalsaveBinsFRVf{1}(edge2plot_vf,1),sumResp/count,'k','LineWidth',2)
box off
set(gcf,'color','w')
xlabel('Vf deg/sec')
ylabel('spikes/sec')

sumResp = zeros(length(totalsaveBinsVmVf{1}(edge2plot_vf,1)),1);
count = 0;
figure(2);clf;
for cells = 1:length(totalsaveBinsVmVf)
    if sum(isnan(totalsaveBinsVmVf{cells}(edge2plot_vf,2))) > 0
    else
        count = count + 1;
        plot(totalsaveBinsVmVf{1}(edge2plot_vf,1),normalize(totalsaveBinsVmVf{cells}(edge2plot_vf,2)),'color',[0.75 0.75 0.75],'LineWidth',1)
        sumResp = sumResp + normalize(totalsaveBinsVmVf{cells}(edge2plot_vf,2));
        hold on
    end
end
plot(totalsaveBinsVmVf{1}(edge2plot_vf,1),sumResp/count,'k','LineWidth',2)
box off
set(gcf,'color','w')
xlabel('Vf deg/sec')
ylabel('Normalized Vm')


sumResp = zeros(length(totalsaveBinsVmVy{1}(edge2plot_vy,1)),1);
count = 0;
figure(3);clf;
for cells = 1:length(totalsaveBinsFRVy)
    if sum(isnan(totalsaveBinsFRVy{cells}(edge2plot_vy,2))) > 0
    else
        count = count + 1;
        plot(totalsaveBinsFRVy{1}(edge2plot_vy,1),normalize(totalsaveBinsFRVy{cells}(edge2plot_vy,2)),'color',[0.75 0.75 0.75],'LineWidth',1)
        sumResp = sumResp + normalize(totalsaveBinsFRVy{cells}(edge2plot_vy,2));
        hold on
    end
end
plot(totalsaveBinsFRVy{1}(edge2plot_vy,1),sumResp/count,'k','LineWidth',2)
box off
set(gcf,'color','w')
xlabel('Vy deg/sec')
ylabel('spikes/sec')

sumResp = zeros(length(totalsaveBinsVmVy{1}(edge2plot_vy,1)),1);
count = 0;
figure(4);clf;
for cells = 1:length(totalsaveBinsVmVy)
    if sum(isnan(totalsaveBinsVmVy{cells}(edge2plot_vy,2))) > 0
    else
        count = count + 1;
        plot(totalsaveBinsVmVy{1}(edge2plot_vy,1),normalize(totalsaveBinsVmVy{cells}(edge2plot_vy,2)),'color',[0.75 0.75 0.75],'LineWidth',1)
        sumResp = sumResp + normalize(totalsaveBinsVmVy{cells}(edge2plot_vy,2));
        hold on
    end
end
plot(totalsaveBinsVmVy{1}(edge2plot_vy,1),sumResp/count,'k','LineWidth',2)
box off
set(gcf,'color','w')
xlabel('Vy deg/sec')
ylabel('Normalized Vm')
    

