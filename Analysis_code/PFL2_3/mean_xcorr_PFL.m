function h = mean_xcorr_PFL(summaryData, t, tData_all, bData_all, savePlots) 
flies = unique(summaryData.Folder); 

aveFlyVmVel = zeros(2001,5,length(flies));
aveFlyFRVel = zeros(2001,5,length(flies));

for f = 1:length(flies)

    flySumData = summaryData(strcmp(summaryData.Folder,flies(f)),:); 
    saveBinsVmVel = zeros(2001,5,size(flySumData,1));
    saveBinsFRVel = zeros(2001,5,size(flySumData,1));

    for chunk = 1:size(flySumData,1)
        tData = tData_all(flySumData.Indices{chunk},:);
        bData = bData_all(flySumData.Indices{chunk},:);
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
        nanIDX  = find( isnan(Vf) | isnan(Vy) | isnan(Vs) );
        % 
        if ~isempty(nanIDX)
            if nanIDX(1) == 1
                Vf = Vf(nanIDX(end)+1:length(Vf)); 
                Vs = Vs(nanIDX(end)+1:length(Vs));
                Vy = Vy(nanIDX(end)+1:length(Vy));
                smVm = smVm(nanIDX(end)+1:length(smVm)); 
                fRate = fRate(nanIDX(end)+1:length(fRate)); 
            end
            
            if nanIDX(end) == length(Vf)
                Vf = Vf(1:nanIDX(1)-1);
                Vs = Vs(1:nanIDX(1)-1);
                Vy = Vy(1:nanIDX(1)-1);
                smVm = smVm(1:nanIDX(1)-1); 
                fRate = fRate(1:nanIDX(1)-1); 
            end
        end

        VyR = Vy;
        Sy = abs(Vy); 
        VsR = Vs;
        Ss = abs(Vs); 


            [cVm, ~] = xcorr(smVm,Vf,1000);
            [cFR, ~] = xcorr(fRate,Vf,1000);

            saveBinsVmVel(:,1,chunk) = (cVm - min(cVm))/(max(cVm) - min(cVm));
            saveBinsFRVel(:,1,chunk) = (cFR - min(cFR))/(max(cFR) - min(cFR));

            [cVm, ~] = xcorr(smVm,VyR,1000);
            [cFR, ~] = xcorr(fRate,VyR,1000);

            saveBinsVmVel(:,2,chunk) = (cVm - min(cVm))/(max(cVm) - min(cVm));
            saveBinsFRVel(:,2,chunk) = (cFR - min(cFR))/(max(cFR) - min(cFR));

            [cVm, ~] = xcorr(smVm,Sy,1000);
            [cFR, ~] = xcorr(fRate,Sy,1000);

            saveBinsVmVel(:,3,chunk) = (cVm - min(cVm))/(max(cVm) - min(cVm));
            saveBinsFRVel(:,3,chunk) = (cFR - min(cFR))/(max(cFR) - min(cFR));

            [cVm, ~] = xcorr(smVm,VsR,1000);
            [cFR, ~] = xcorr(fRate,VsR,1000);

            saveBinsVmVel(:,4,chunk) = (cVm - min(cVm))/(max(cVm) - min(cVm));
            saveBinsFRVel(:,4,chunk) = (cFR - min(cFR))/(max(cFR) - min(cFR));

            [cVm, lagsVm] = xcorr(smVm,Ss,1000);
            [cFR, lagsFR] = xcorr(fRate,Ss,1000);

            saveBinsVmVel(:,5,chunk) = (cVm - min(cVm))/(max(cVm) - min(cVm));
            saveBinsFRVel(:,5,chunk) = (cFR - min(cFR))/(max(cFR) - min(cFR));
    end
    
    aveFlyVmVel(:,:,f) = mean(saveBinsVmVel,3);
    aveFlyFRVel(:,:,f) = mean(saveBinsFRVel,3);
end

if length(flies) ~= 1
    finalAveVm = mean(aveFlyVmVel(f,:,:),1);
    finalAveFR = mean(aveFlyFRVel(f,:,:),1);
else
    finalAveVm = aveFlyVmVel;
    finalAveFR = aveFlyFRVel(:,:,f);
end

h = figure();
set(gcf,'color','w')
set(gcf,'renderer','painters')

if length(flies) == 1
    for c = 1:size(flySumData,1)
        subplot(2,5,1);
        plot(lagsFR,saveBinsFRVel(:,1,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('Vf FR')
    subplot(2,5,2);
        plot(lagsFR,saveBinsFRVel(:,2,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('Vy FR')  
    subplot(2,5,3);
        plot(lagsFR,saveBinsFRVel(:,3,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('yaw speed FR')

    subplot(2,5,4);
        plot(lagsFR,saveBinsFRVel(:,4,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('Vs FR')
    subplot(2,5,5);
        plot(lagsFR,saveBinsFRVel(:,5,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('side speed FR')
        %%
    subplot(2,5,6);
        plot(lagsVm,saveBinsVmVel(:,1,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('Vf Vm')
    subplot(2,5,7);
        plot(lagsVm,saveBinsVmVel(:,2,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('Vy Vm')  
    subplot(2,5,8);
        plot(lagsVm,saveBinsVmVel(:,3,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('yaw speed Vm')
    subplot(2,5,9);
        plot(lagsVm,saveBinsVmVel(:,4,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('Vs Vm')
    subplot(2,5,10);
        plot(lagsVm,saveBinsVmVel(:,5,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('side speed Vm')
    end
else
    for c = 1:length(flies)
    subplot(2,5,1);
        plot(lagsFR,aveFlyFRVel(:,1,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('Vf FR')
    subplot(2,5,2);
        plot(lagsFR,aveFlyFRVel(:,2,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('Vy FR')  
    subplot(2,5,3);
        plot(lagsFR,aveFlyFRVel(:,3,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('yaw speed FR')

    subplot(2,5,4);
        plot(lagsFR,aveFlyFRVel(:,4,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('Vs FR')
    subplot(2,5,5);
        plot(lagsFR,aveFlyFRVel(:,5,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('side speed FR')
        %%
    subplot(2,5,6);
        plot(lagsVm,aveFlyVmVel(:,1,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('Vf Vm')
    subplot(2,5,7);
        plot(lagsVm,aveFlyVmVel(:,2,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('Vy Vm')  
    subplot(2,5,8);
        plot(lagsVm,aveFlyVmVel(:,3,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('yaw speed Vm')
    subplot(2,5,9);
        plot(lagsVm,aveFlyVmVel(:,4,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('Vs Vm')
    subplot(2,5,10);
        plot(lagsVm,aveFlyVmVel(:,5,c),'color',[0.75,0.75,0.75])
        hold on
        legend('hide')
        title('side speed Vm')
    end
end

 subplot(2,5,1);
        plot(lagsFR,finalAveFR(:,1),'k','LineWidth',1.5)
        hold on
        legend('hide')
        title('Vf FR')
    subplot(2,5,2);
        plot(lagsFR,finalAveFR(:,2),'k','LineWidth',1.5)
        hold on
        legend('hide')
        title('Vy FR')  
    subplot(2,5,3);
        plot(lagsFR,finalAveFR(:,3),'k','LineWidth',1.5)
        hold on
        legend('hide')
        title('yaw speed FR')

    subplot(2,5,4);
        plot(lagsFR,finalAveFR(:,4),'k','LineWidth',1.5)
        hold on
        legend('hide')
        title('Vs FR')
    subplot(2,5,5);
        plot(lagsFR,finalAveFR(:,5),'k','LineWidth',1.5)
        hold on
        legend('hide')
        title('side speed FR')
        %%
    subplot(2,5,6);
        plot(lagsVm,finalAveVm(:,1),'k','LineWidth',1.5)
        hold on
        legend('hide')
        title('Vf Vm')
    subplot(2,5,7);
        plot(lagsVm,finalAveVm(:,2),'k','LineWidth',1.5)
        hold on
        legend('hide')
        title('Vy Vm')  
    subplot(2,5,8);
        plot(lagsVm,finalAveVm(:,3),'k','LineWidth',1.5)
        hold on
        legend('hide')
        title('yaw speed Vm')
    subplot(2,5,9);
        plot(lagsVm,finalAveVm(:,4),'k','LineWidth',1.5)
        hold on
        legend('hide')
        title('Vs Vm')
    subplot(2,5,10);
        plot(lagsVm,finalAveVm(:,5),'k','LineWidth',1.5)
        hold on
        legend('hide')
        title('side speed Vm')


%%

    if savePlots
        saveas(gcf,fullfile(folder, 'figures',['trial_',num2str(t),'_average_xcorr.fig']))
    end