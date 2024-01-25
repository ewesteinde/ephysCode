function h = basic_xcorr_PFL(folder, t, tData, bData, savePlots) 

saveBinsVmVel = cell(1,5);
saveBinsFRVel = cell(1,5); 



%seemingly randomly some exps have a continuous stream of NaN values at
%start or end of trials, will need to cut these sections in order to run

Vf = bData.vel_for;
Vs = bData.vel_side;
Vy = bData.vel_yaw;
try
    smVm = tData.smooth_Vm;
catch
    smVm = tData.smoothVm;
end
fRate = tData.fRate_sec;
nanIDX  = find( isnan(Vf) );
% 
if ~isempty(nanIDX)
    if nanIDX(1) == 1
        Vf = Vf(nanIDX(end)+1:length(Vf)); 
        Vs = Vs(nanIDX(end)+1:length(Vs));
        Vy = Vy(nanIDX(end)+1:length(Vy));
        smVm = smVm(nanIDX(end)+1:length(smVm)); 
        fRate = fRate(nanIDX(end)+1:length(fRate));  
    elseif nanIDX(end) == length(Vf)
        Vf = Vf(1:nanIDX(1)-1);
        Vs = Vs(1:nanIDX(1)-1);
        Vy = Vy(1:nanIDX(1)-1);
        smVm = smVm(1:nanIDX(1)-1); 
        fRate = fRate(1:nanIDX(1)-1); 
    else
        error('NaN values inside velDot')
    end
end

nanIDX  = find( isnan(Vy) );
% 
if ~isempty(nanIDX)
    if nanIDX(1) == 1
        Vf = Vf(nanIDX(end)+1:length(Vf)); 
        Vs = Vs(nanIDX(end)+1:length(Vs));
        Vy = Vy(nanIDX(end)+1:length(Vy));
        smVm = smVm(nanIDX(end)+1:length(smVm)); 
        fRate = fRate(nanIDX(end)+1:length(fRate));  
    elseif nanIDX(end) == length(Vf)
        Vf = Vf(1:nanIDX(1)-1);
        Vs = Vs(1:nanIDX(1)-1);
        Vy = Vy(1:nanIDX(1)-1);
        smVm = smVm(1:nanIDX(1)-1); 
        fRate = fRate(1:nanIDX(1)-1); 
    else
        error('NaN values inside velDot')
    end
end

VyR = Vy;
Sy = abs(Vy); 
VsR = Vs;
Ss = abs(Vs); 

  
    [cVm, lagsVm] = xcorr(smVm,Vf,1000);
    [cFR, lagsFR] = xcorr(fRate,Vf,1000);

    saveBinsVmVel{1} = [(lagsVm./1000)', cVm];
    saveBinsFRVel{1} = [(lagsFR./1000)', cFR];
    
    [cVm, lagsVm] = xcorr(smVm,VyR,1000);
    [cFR, lagsFR] = xcorr(fRate,VyR,1000);

    saveBinsVmVel{2} = [(lagsVm./1000)', cVm];
    saveBinsFRVel{2} = [(lagsFR./1000)', cFR];
    
    [cVm, lagsVm] = xcorr(smVm,Sy,1000);
    [cFR, lagsFR] = xcorr(fRate,Sy,1000);

    saveBinsVmVel{3} = [(lagsVm./1000)', cVm];
    saveBinsFRVel{3} = [(lagsFR./1000)', cFR];
    
    [cVm, lagsVm] = xcorr(smVm,VsR,1000);
    [cFR, lagsFR] = xcorr(fRate,VsR,1000);

    saveBinsVmVel{4} = [(lagsVm./1000)', cVm];
    saveBinsFRVel{4} = [(lagsFR./1000)', cFR];
    
    [cVm, lagsVm] = xcorr(smVm,Ss,1000);
    [cFR, lagsFR] = xcorr(fRate,Ss,1000);

    saveBinsVmVel{5} = [(lagsVm./1000)', cVm];
    saveBinsFRVel{5} = [(lagsFR./1000)', cFR];

h= figure(); clf; 

subplot(2,5,1);
    plot(saveBinsFRVel{1}(:,1),saveBinsFRVel{1}(:,2))
    xline(saveBinsFRVel{1}(saveBinsFRVel{1}(:,2) == max(saveBinsFRVel{1}(:,2)),1),'-',num2str(saveBinsFRVel{1}(saveBinsFRVel{1}(:,2) == max(saveBinsFRVel{1}(:,2)),1)))
    xline(saveBinsFRVel{1}(saveBinsFRVel{1}(:,2) == min(saveBinsFRVel{1}(:,2)),1),'-',num2str(saveBinsFRVel{1}(saveBinsFRVel{1}(:,2) == min(saveBinsFRVel{1}(:,2)),1)))
    legend('hide')
    title('Vf FR')
subplot(2,5,2);
    plot(saveBinsFRVel{2}(:,1),saveBinsFRVel{2}(:,2))
    xline(saveBinsFRVel{2}(saveBinsFRVel{2}(:,2) == max(saveBinsFRVel{2}(:,2)),1),'-',num2str(saveBinsFRVel{2}(saveBinsFRVel{2}(:,2) == max(saveBinsFRVel{2}(:,2)),1)))
    xline(saveBinsFRVel{2}(saveBinsFRVel{2}(:,2) == min(saveBinsFRVel{2}(:,2)),1),'-',num2str(saveBinsFRVel{2}(saveBinsFRVel{2}(:,2) == min(saveBinsFRVel{2}(:,2)),1)))
    legend('hide')
    title('Vy FR')  
subplot(2,5,3);
    plot(saveBinsFRVel{3}(:,1),saveBinsFRVel{3}(:,2))
    xline(saveBinsFRVel{3}(saveBinsFRVel{3}(:,2) == max(saveBinsFRVel{3}(:,2)),1),'-',num2str(saveBinsFRVel{3}(saveBinsFRVel{3}(:,2) == max(saveBinsFRVel{3}(:,2)),1)))
    legend('hide')
    title('yaw speed FR')

subplot(2,5,4);
    plot(saveBinsFRVel{4}(:,1),saveBinsFRVel{4}(:,2))
    xline(saveBinsFRVel{4}(saveBinsFRVel{4}(:,2) == max(saveBinsFRVel{4}(:,2)),1),'-',num2str(saveBinsFRVel{4}(saveBinsFRVel{4}(:,2) == max(saveBinsFRVel{4}(:,2)),1)))
    xline(saveBinsFRVel{4}(saveBinsFRVel{4}(:,2) == min(saveBinsFRVel{4}(:,2)),1),'-',num2str(saveBinsFRVel{4}(saveBinsFRVel{4}(:,2) == min(saveBinsFRVel{4}(:,2)),1)))
    legend('hide')
    title('Vs FR')
subplot(2,5,5);
    plot(saveBinsFRVel{5}(:,1),saveBinsFRVel{5}(:,2))
    xline(saveBinsFRVel{5}(saveBinsFRVel{5}(:,2) == max(saveBinsFRVel{5}(:,2)),1),'-',num2str(saveBinsFRVel{5}(saveBinsFRVel{5}(:,2) == max(saveBinsFRVel{5}(:,2)),1)))
    legend('hide')
    title('side speed FR')
    %%
subplot(2,5,6);
    plot(saveBinsVmVel{1}(:,1),saveBinsVmVel{1}(:,2))
    xline(saveBinsVmVel{1}(saveBinsVmVel{1}(:,2) == max(saveBinsVmVel{1}(:,2)),1),'-',num2str(saveBinsVmVel{1}(saveBinsVmVel{1}(:,2) == max(saveBinsVmVel{1}(:,2)),1)))
    xline(saveBinsVmVel{1}(saveBinsVmVel{1}(:,2) == min(saveBinsVmVel{1}(:,2)),1),'-',num2str(saveBinsVmVel{1}(saveBinsVmVel{1}(:,2) == min(saveBinsVmVel{1}(:,2)),1)))
    legend('hide')
    title('Vf Vm')
subplot(2,5,7);
    plot(saveBinsVmVel{2}(:,1),saveBinsVmVel{2}(:,2))
    xline(saveBinsVmVel{2}(saveBinsVmVel{2}(:,2) == max(saveBinsVmVel{2}(:,2)),1),'-',num2str(saveBinsVmVel{2}(saveBinsVmVel{2}(:,2) == max(saveBinsVmVel{2}(:,2)),1)))
    xline(saveBinsVmVel{2}(saveBinsVmVel{2}(:,2) == min(saveBinsVmVel{2}(:,2)),1),'-',num2str(saveBinsVmVel{2}(saveBinsVmVel{2}(:,2) == min(saveBinsVmVel{2}(:,2)),1)))
    legend('hide')
    title('Vy Vm')  
subplot(2,5,8);
    plot(saveBinsVmVel{3}(:,1),saveBinsVmVel{3}(:,2))
    xline(saveBinsVmVel{3}(saveBinsVmVel{3}(:,2) == max(saveBinsVmVel{3}(:,2)),1),'-',num2str(saveBinsVmVel{3}(saveBinsVmVel{3}(:,2) == max(saveBinsVmVel{3}(:,2)),1)))
    legend('hide')
    title('yaw speed Vm')

subplot(2,5,9);
    plot(saveBinsVmVel{4}(:,1),saveBinsVmVel{4}(:,2))
    xline(saveBinsVmVel{4}(saveBinsVmVel{4}(:,2) == max(saveBinsVmVel{4}(:,2)),1),'-',num2str(saveBinsVmVel{4}(saveBinsVmVel{4}(:,2) == max(saveBinsVmVel{4}(:,2)),1)))
    xline(saveBinsVmVel{4}(saveBinsVmVel{4}(:,2) == min(saveBinsVmVel{4}(:,2)),1),'-',num2str(saveBinsVmVel{4}(saveBinsVmVel{4}(:,2) == min(saveBinsVmVel{4}(:,2)),1)))
    legend('hide')
    title('Vs Vm')
subplot(2,5,10);
    plot(saveBinsVmVel{5}(:,1),saveBinsVmVel{5}(:,2))
    xline(saveBinsVmVel{5}(saveBinsVmVel{5}(:,2) == max(saveBinsVmVel{5}(:,2)),1),'-',num2str(saveBinsVmVel{5}(saveBinsVmVel{5}(:,2) == max(saveBinsVmVel{5}(:,2)),1)))
    legend('hide')
    title('side speed Vm')


%%

    if savePlots
        saveas(gcf,fullfile(folder, 'figures',['trial_',num2str(t),'_wholeTrial_xcorr.fig']))
    end