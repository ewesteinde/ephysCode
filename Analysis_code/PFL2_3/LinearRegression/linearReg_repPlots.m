%% separate idx at which fly is moving in CL vs OL
    
    startOL = 178.189*1000; % sec
    endOL = 203.249*1000;
    lengthCL = 179.99*1000;
    
% start = 1819.45;
% finish = 1844.47;
%     startOL  = find(processed_behaviourData.time == start); 
%     endOL = find(processed_behaviourData.time == finish); 

[CL, OL, CL_startStopIdx, OL_startStopIdx] = separateCL_OL(startOL, endOL, lengthCL, processed_behaviourData, processed_trialData);

%% use a sliding 1 min window w/ overlap
CL_chunks = [];
window = 60; 
window = window * 1000;
step = 30 * 1000;
count = 1; 
for chunk = 1:length(CL_startStopIdx)
    chunk_idx = CL_startStopIdx(chunk, 1):CL_startStopIdx(chunk, 2);
    for i = 1:step:(length(chunk_idx) - window + 1)
        CL_chunks(count,1) = chunk_idx(i);
        CL_chunks(count,2) = chunk_idx(i + window - 1); 
        count = count + 1;
    end
    
    if (i + 2*window - step) > length(chunk_idx)
        CL_chunks(count,1) = chunk_idx(i + window - step);
        CL_chunks(count,2) = chunk_idx(end); 
        count = count + 1;
    end
end

%%
seg = [1:length(CL_chunks)];
[BVm_nobasis, ave_deg_nobasis, ave_rho_nobasis, prefHeads_nobasis, R2_Vm_nobasis] = linearReg_chunks(seg, processed_behaviourData, processed_trialData, CL_chunks, fileName);

%%
chunks = (1:length(BVm(1,:)))/2;
fig = figure(1);clf;
subplot(7,1,1)
    plot(chunks,BVm(1,:),'k')
    ylabel('B Heading')
    box off
subplot(7,1,2);
    plot(chunks,BVm(2,:),'k')
    ylabel('B Forward Vel')
    box off
subplot(7,1,3);
    plot(chunks,BVm(3,:),'k')
    ylabel('B Yaw Vel') 
    box off
subplot(7,1,4);
    plot(chunks,ave_deg,'k')
    ylabel('Heading')
    ylim([-200 200])
    box off
subplot(7,1,5);
    plot(chunks,ave_rho,'k')
    ylabel('Vector Strength')
    ylim([0 1])
    box off
subplot(7,1,6);
    plot(chunks,wrapTo180(prefHeads),'k')
    ylim([-200 200])
    ylabel('Preferred Heading')
    box off
subplot(7,1,7)
    plot(chunks,R2_Vm,'k')
    ylabel('R^2')
    xlabel('time (min)')
set(gcf,'color','white')
box off 

%%

fig = figure(2);clf;
left_color = [0 0 0];
right_color = [1 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
    plot(chunks,R2_Vm_basis,'k')
    ylabel('R^2 heading & vf')
    box off
yyaxis right
    plot(chunks,R2_Vm,'r')
    ylabel('R^2 only heading')
    xlabel('Time (min)')
     ylim([0 0.5])
    box off
%      y = min(prefHeads);
%      l = line([chunks(end)-3 chunks(end)-1], [y y]);
%      l.Color = 'k';
% set(gca,'XColor','none')
set(gcf,'color',[1 1 1])