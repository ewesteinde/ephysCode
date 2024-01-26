function [vf_coeff, vy_coeff] = PFL2_3_behaviourVSactivity_lineplots_new(prefHead, step, range, processed_trialData, processed_behaviourData) 

    angles = {num2str(wrapTo180(prefHead+180)), strcat(num2str(wrapTo180(prefHead+90)),'/',num2str(wrapTo180(prefHead-90))), num2str(prefHead)};
    headings = [wrapTo180(prefHead+180),wrapTo180(prefHead+90), prefHead];
    

        fRate_sec = processed_trialData.fRate_sec;
        try
            Vm = processed_trialData.smooth_Vm;
        catch
            Vm = processed_trialData.smoothVm;
        end
        angle = processed_behaviourData.angle;
        Vf = processed_behaviourData.vel_for;
        Vs = processed_behaviourData.vel_side;
        Vy = processed_behaviourData.vel_yaw;

    
    speed = sqrt(Vs.^2 + Vf.^2); 

    count = 1; 
    saveBinsVmS = cell(1,3);
    saveBinsFRS = cell(1,3);
    saveBinsVmVf = cell(1,3);
    saveBinsFRVf = cell(1,3); 
    saveBinsVmVs = cell(1,3);
    saveBinsFRVs = cell(1,3);
    saveBinsVmVy = cell(1,3);
    saveBinsFRVy = cell(1,3);
    
    for head = headings
        VfinHead = []; 
        VminHead = []; 
        fRateinHead = []; 
        index = [];

        if count == 2
            rangeh = range/2;
            
            lim1 = head-rangeh/2;
            lim2 = head+rangeh/2;
            
            lim3 = wrapTo180(head+180)-rangeh/2;
            lim4 = wrapTo180(head+180)+rangeh/2;
        
            if abs(lim1) > 180
                [index1] = find(angle >= wrapTo180(lim1) | angle <= lim2);
                upperlim1 = lim2;
                lowerlim1 = wrapTo180(lim1);
            elseif abs(lim2) > 180
                [index1] = find(angle >= lim1 | angle <= wrapTo180(lim2));
                upperlim1 = wrapTo180(lim2);
                lowerlim1 = lim1;
            else
                [index1] = find(angle >= lim1 & angle <= lim2);
                upperlim1 = lim2;
                lowerlim1 = lim1;
            end
            
            if abs(lim3) > 180
                [index2] = find(angle >= wrapTo180(lim3) | angle <= lim4);
                upperlim2 = lim4;
                lowerlim2 = wrapTo180(lim3);
            elseif abs(lim4) > 180
                [index2] = find(angle >= lim3 | angle <= wrapTo180(lim4));
                upperlim2 = wrapTo180(lim4);
                lowerlim2 = lim3;
            else
                [index2] = find(angle >= lim3 & angle <= lim4);
                upperlim2 = lim4;
                lowerlim2 = lim3;
            end
            
            index = cat(1,index1, index2); 
        else

            lim1 = head-range/2;
            lim2 = head+range/2;
        
            if abs(lim1) > 180
                [index] = find(angle >= wrapTo180(lim1) | angle <= lim2);
                upperlim = lim2;
                lowerlim = wrapTo180(lim1);
            elseif abs(lim2) > 180
                [index] = find(angle >= lim1 | angle <= wrapTo180(lim2));
                upperlim = wrapTo180(lim2);
                lowerlim = lim1;
            else
                [index] = find(angle >= lim1 & angle <= lim2);
                upperlim = lim2;
                lowerlim = lim1;
            end
        end
        
        VfinHead = Vf(index); 
        VsinHead = Vs(index);
        VyinHead = Vy(index);
        SinHead = speed(index);
        VminHead = Vm(index); 
        fRateinHead = fRate_sec(index);
        
        no0vel_idx = find(SinHead > 1.5);
        
        VfinHead = Vf(no0vel_idx); 
        VsinHead = Vs(no0vel_idx);
        VyinHead = Vy(no0vel_idx);
        SinHead = speed(no0vel_idx);
        VminHead = Vm(no0vel_idx); 
        fRateinHead = fRate_sec(no0vel_idx);

%                 edges = [-(step/2):step:max(SinHead(:,1))]; %start at -step/2 so center of first bin is 0mm/s
%                 [N, edges, bin] = histcounts(SinHead(:,1), edges);
%                 tempFR = accumarray(bin+1, fRateinHead, [length(edges) 1]);
%                 tempVm = accumarray(bin+1, VminHead, [length(edges) 1]);
%                 mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
%                 mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
%                 centers = edges(1:end-1)+diff(edges)/2;
%                 
%                 saveBinsVmS{count} = [centers' mean_binVm];
%                 saveBinsFRS{count} = [centers' mean_binFR];
                
                
                edges = [-4:0.5:8]; 
                [N, edges, bin] = histcounts(VfinHead(:,1), edges);
                tempFR = accumarray(bin+1, fRateinHead, [length(edges) 1]);
                tempVm = accumarray(bin+1, VminHead, [length(edges) 1]);
                mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
                mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
                centers = edges(1:end-1)+diff(edges)/2;
                
                saveBinsVmVf{count} = [centers' mean_binVm];
                saveBinsFRVf{count} = [centers' mean_binFR];
                
                edges = [-200:10:200]; 
                [N, edges, bin] = histcounts(VyinHead(:,1), edges);
                tempFR = accumarray(bin+1, fRateinHead, [length(edges) 1]);
                tempVm = accumarray(bin+1, VminHead, [length(edges) 1]);
                mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
                mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
                centers = edges(1:end-1)+diff(edges)/2;
                
                saveBinsVmVy{count} = [centers' mean_binVm];
                saveBinsFRVy{count} = [centers' mean_binFR];
%                 
%                 edges = [min(VsinHead(:,1)):step:max(VsinHead(:,1))]; 
%                 [N, edges, bin] = histcounts(VsinHead(:,1), edges);
%                 tempFR = accumarray(bin+1, fRateinHead, [length(edges) 1]);
%                 tempVm = accumarray(bin+1, VminHead, [length(edges) 1]);
%                 mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
%                 mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
%                 centers = edges(1:end-1)+diff(edges)/2;
%                 
%                 saveBinsVmVs{count} = [centers' mean_binVm];
%                 saveBinsFRVs{count} = [centers' mean_binFR];
        count = count + 1; 
    end
    
%%

g = figure(); clf; 

   subplot(2,1,1);
   hold on
    plot(saveBinsVmVf{1}(:,1),saveBinsVmVf{1}(:,2), 'color', [0.5 0.5 0.5])  
    vf_antiPrefHead_coeff = polyfit(saveBinsVmVf{1}(~isnan(saveBinsVmVf{1}(:,2)),1), saveBinsVmVf{1}(~isnan(saveBinsVmVf{1}(:,2)),2), 1);
    xFit = linspace(min(saveBinsVmVf{1}(:,1)), max(saveBinsVmVf{1}(:,1)), 1000);
    yFit = polyval(vf_antiPrefHead_coeff , xFit);
    plot(xFit, yFit,'color',[0.5 0.5 0.5])
    plot(saveBinsVmVf{2}(:,1),saveBinsVmVf{2}(:,2),'color', [0.5 0 0])
    plot(saveBinsVmVf{3}(:,1),saveBinsVmVf{3}(:,2),'color', [1 0 0])
    vf_prefHead_coeff = polyfit(saveBinsVmVf{3}(~isnan(saveBinsVmVf{3}(:,2)),1), saveBinsVmVf{3}(~isnan(saveBinsVmVf{3}(:,2)),2), 1);
    xFit = linspace(min(saveBinsVmVf{3}(:,1)), max(saveBinsVmVf{3}(:,1)), 1000);
    yFit = polyval(vf_prefHead_coeff , xFit);
    plot(xFit, yFit,'r')

   % xlim([minValPlot+1 maxValPlot-1])
    xlabel('Vf mm/sec')
    ylabel('Vm (mV)')
    xlim([-2 10])
    set(gcf,'color',[1 1 1])
    box off
    
%     subplot(2,1,2);
%     plot(saveBinsFRVf{1}(:,1),saveBinsFRVf{1}(:,2), 'color', [0.5 0.5 0.5])
%     hold on
%     plot(saveBinsFRVf{2}(:,1),saveBinsFRVf{2}(:,2),'color', [0.5 0 0])
%     plot(saveBinsFRVf{3}(:,1),saveBinsFRVf{3}(:,2),'color', [1 0 0])
% 
%    % xlim([minValPlot+1 maxValPlot-1])
%     xlabel('Vf mm/sec')
%     ylabel('spikes/sec')
%     xlim([-2 10])
%     set(gcf,'color',[1 1 1])
%     box off
    
% yaw velocity 
    subplot(2,1,2);
    hold on
    plot(saveBinsVmVy{1}(:,1),saveBinsVmVy{1}(:,2), 'color', [0.5 0.5 0.5])  
    vy_antiPrefHead_coeff = polyfit(saveBinsVmVy{1}(~isnan(saveBinsVmVy{1}(:,2)),1), saveBinsVmVy{1}(~isnan(saveBinsVmVy{1}(:,2)),2), 1);
    xFit = linspace(min(saveBinsVmVy{1}(:,1)), max(saveBinsVmVy{1}(:,1)), 1000);
    yFit = polyval(vy_antiPrefHead_coeff , xFit);
    plot(xFit, yFit,'color',[0.5 0.5 0.5])
    plot(saveBinsVmVy{2}(:,1),saveBinsVmVy{2}(:,2),'color', [0.5 0 0])
    plot(saveBinsVmVy{3}(:,1),saveBinsVmVy{3}(:,2),'color', [1 0 0])
    vy_prefHead_coeff = polyfit(saveBinsVmVy{3}(~isnan(saveBinsVmVy{3}(:,2)),1), saveBinsVmVy{3}(~isnan(saveBinsVmVy{3}(:,2)),2), 1);
    xFit = linspace(min(saveBinsVmVy{3}(:,1)), max(saveBinsVmVy{3}(:,1)), 1000);
    yFit = polyval(vy_prefHead_coeff , xFit);
    plot(xFit, yFit,'r')

    %xlim([minValPlot+1 maxValPlot-1])
    xlabel('Vy mm/sec')
    ylabel('Vm (mV)')
    xlim([-250 250])
    set(gcf,'color',[1 1 1])
    box off
    
%     subplot(2,1,2);
%     plot(saveBinsFRVy{1}(:,1),saveBinsFRVy{1}(:,2), 'color', [0.5 0.5 0.5])
%     hold on
%     plot(saveBinsFRVy{2}(:,1),saveBinsFRVy{2}(:,2),'color', [0.5 0 0])
%     plot(saveBinsFRVy{3}(:,1),saveBinsFRVy{3}(:,2),'color', [1 0 0])
% 
%     %xlim([minValPlot+1 maxValPlot-1])
%     xlabel('Vy mm/sec')
%     ylabel('spikes/sec')
%         xlim([-250 250])
%     set(gcf,'color',[1 1 1])
%     box off
    
vf_coeff(1,:) = vf_prefHead_coeff;
vf_coeff(2,:) = vf_antiPrefHead_coeff;

vy_coeff(1,:) = vy_prefHead_coeff;
vy_coeff(2,:) = vy_antiPrefHead_coeff;
    


end
