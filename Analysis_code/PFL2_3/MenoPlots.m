figure()
plot(processed_behaviourData.angle)

meno_idx = 239992:length(CL.angle); % 265070:length(processed_behaviourData.angle);%
nMeno_idx = 1:239992;




figure();clf;
Vm = CL.smooth_Vm(meno_idx);
angle = CL.angle(meno_idx);
edges = [-180:30:180];
[N,centers, mean_bin] = create_binned_mean(angle, Vm, edges);
% centers(end + 1) = centers(1);
% mean_bin(end+1) = mean_bin(1);
plot(centers, mean_bin,'k');
% polarplot(wrapTo2Pi(deg2rad(centers)),mean_bin,'k')
% rlim([-61,-51])
xlabel('Cue Position');
ylabel('Vm'); 
title('Meno heading pref')
set(gcf,'Renderer', 'painters','color','w');
box off
xlim([-180,180])
ylim([-61,-51])

%[pval_meno,z_meno] = circ_rtest(wrapTo2Pi(deg2rad(centers)),N,15);

figure();clf;
Vm = CL.smooth_Vm(nMeno_idx);
angle = CL.angle(nMeno_idx);
edges = [-180:30:180];
[N,centers, mean_bin] = create_binned_mean(angle, Vm, edges);
% centers(end + 1) = centers(1);
% mean_bin(end+1) = mean_bin(1);
%polarplot(wrapTo2Pi(deg2rad(centers)),mean_bin,'k')
%rlim([-61,-51])
plot(centers, mean_bin,'k');
xlabel('Cue Position');
ylabel('Vm'); 
title('nMeno heading pref')
set(gcf,'Renderer', 'painters','color','w');
grid off
%[pval_nMeno,z_nMeno] = circ_rtest(mean_bin,N,15);
box off
xlim([-180,180])
ylim([-61,-51])

figure();clf;
Vm = processed_trialData.smooth_Vm(meno_idx);
angle = processed_behaviourData.angle(meno_idx);
edges = [-180:20:180];
[N,centers, mean_bin] = create_binned_mean(angle, Vm, edges);
% centers(end + 1) = centers(1);
% mean_bin(end+1) = mean_bin(1);
plot(centers, mean_bin,'k');
% polarplot(wrapTo2Pi(deg2rad(centers)),mean_bin,'k')
% rlim([-61,-51])
xlabel('Cue Position');
ylabel('Vm'); 
title('Meno heading pref')
set(gcf,'Renderer', 'painters','color','w');
box off
xlim([-180,180])
ylim([-61,-51])

%[pval_meno,z_meno] = circ_rtest(wrapTo2Pi(deg2rad(centers)),N,15);

figure();clf;
Vm = processed_trialData.smooth_Vm(nMeno_idx);
angle = processed_behaviourData.angle(nMeno_idx);
edges = [-180:20:180];
[N,centers, mean_bin] = create_binned_mean(angle, Vm, edges);
% centers(end + 1) = centers(1);
% mean_bin(end+1) = mean_bin(1);
%polarplot(wrapTo2Pi(deg2rad(centers)),mean_bin,'k')
%rlim([-61,-51])
plot(centers, mean_bin,'k');
xlabel('Cue Position');
ylabel('Vm'); 
title('nMeno heading pref')
set(gcf,'Renderer', 'painters','color','w');
grid off
%[pval_nMeno,z_nMeno] = circ_rtest(mean_bin,N,15);
box off
xlim([-180,180])
ylim([-61,-51])