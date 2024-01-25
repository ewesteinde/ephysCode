function Ionto_behaviourVSactivity_lineplots(folder, prefHead, range, processed_trialData, processed_behaviourData, savePlot) 

%% 
    %name = strcat(date,', ',cell_num,', ',trial,' prefAngle: ',num2str(prefHead),' BinSize: ',num2str(range));
    angles = {num2str(wrapTo180(prefHead+180)), strcat(num2str(wrapTo180(prefHead+90)),'/',num2str(wrapTo180(prefHead-90))), num2str(prefHead)};
    headings = [wrapTo180(prefHead+180),wrapTo180(prefHead+90), prefHead];
    
    total_mov_mm = abs(processed_behaviourData.vel_for + abs(processed_behaviourData.vel_side) + abs(deg2rad(processed_behaviourData.vel_yaw))*4.5);
    no0vel_idx = find(total_mov_mm > 3);

        fRate_sec = processed_trialData.fRate_sec(no0vel_idx);
        Vm = processed_trialData.smoothVm(no0vel_idx);
        angle = processed_behaviourData.angle(no0vel_idx);
        Vf = processed_behaviourData.vel_for(no0vel_idx);
        Vs = abs(processed_behaviourData.vel_side(no0vel_idx));
        Vy = abs(processed_behaviourData.vel_yaw(no0vel_idx));

    
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

                edges = [-(1/2):1:max(SinHead(:,1))]; %start at -step/2 so center of first bin is 0mm/s
                [N, edges, bin] = histcounts(SinHead(:,1), edges);
                tempFR = accumarray(bin+1, fRateinHead, [length(edges) 1]);
                tempVm = accumarray(bin+1, VminHead, [length(edges) 1]);
                mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
                mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
                centers = edges(1:end-1)+diff(edges)/2;
                
                saveBinsVmS{count} = [centers' mean_binVm];
                saveBinsFRS{count} = [centers' mean_binFR];
                
                
                edges = [-5:1:10]; 
                [N, edges, bin] = histcounts(VfinHead(:,1), edges);
                tempFR = accumarray(bin+1, fRateinHead, [length(edges) 1]);
                tempVm = accumarray(bin+1, VminHead, [length(edges) 1]);
                mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
                mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
                centers = edges(1:end-1)+diff(edges)/2;
                
                saveBinsVmVf{count} = [centers' mean_binVm];
                saveBinsFRVf{count} = [centers' mean_binFR];
                
                edges = [0:5:200]; 
                [N, edges, bin] = histcounts(VyinHead(:,1), edges);
                tempFR = accumarray(bin+1, fRateinHead, [length(edges) 1]);
                tempVm = accumarray(bin+1, VminHead, [length(edges) 1]);
                mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
                mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
                centers = edges(1:end-1)+diff(edges)/2;
                
                saveBinsVmVy{count} = [centers' mean_binVm];
                saveBinsFRVy{count} = [centers' mean_binFR];
                
                edges = [0:0.5:8]; 
                [N, edges, bin] = histcounts(VsinHead(:,1), edges);
                tempFR = accumarray(bin+1, fRateinHead, [length(edges) 1]);
                tempVm = accumarray(bin+1, VminHead, [length(edges) 1]);
                mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
                mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
                centers = edges(1:end-1)+diff(edges)/2;
                
                saveBinsVmVs{count} = [centers' mean_binVm];
                saveBinsFRVs{count} = [centers' mean_binFR];
        count = count + 1; 
    end
    
%%

g = figure(); clf; 

   subplot(2,1,1);
    plot(saveBinsVmVf{1}(:,1),saveBinsVmVf{1}(:,2), 'color', [0.5 0.5 0.5])  
    hold on 
    plot(saveBinsVmVf{2}(:,1),saveBinsVmVf{2}(:,2),'color', [0.5 0 0])
    plot(saveBinsVmVf{3}(:,1),saveBinsVmVf{3}(:,2),'color', [1 0 0])

   % xlim([minValPlot+1 maxValPlot-1])
    xlabel('Vf mm/sec')
    ylabel('Vm (mV)')
    set(gcf,'color',[1 1 1])
    box off
    
    subplot(2,1,2);
    plot(saveBinsFRVf{1}(:,1),saveBinsFRVf{1}(:,2), 'color', [0.5 0.5 0.5])
    hold on
    plot(saveBinsFRVf{2}(:,1),saveBinsFRVf{2}(:,2),'color', [0.5 0 0])
    plot(saveBinsFRVf{3}(:,1),saveBinsFRVf{3}(:,2),'color', [1 0 0])

   % xlim([minValPlot+1 maxValPlot-1])
    xlabel('Vf mm/sec')
    ylabel('spikes/sec')
    set(gcf,'color',[1 1 1])
    box off
    
% sideways velocity 

%     subplot(4,2,2);
%     plot(saveBinsVmVs{1}(:,1),saveBinsVmVs{1}(:,2), 'color', [0.5 0.5 0.5])  
%     hold on 
%     plot(saveBinsVmVs{2}(:,1),saveBinsVmVs{2}(:,2),'color', [0.5 0 0])
%     plot(saveBinsVmVs{3}(:,1),saveBinsVmVs{3}(:,2),'color', [1 0 0])
% 
%    % xlim([minValPlot+1 maxValPlot-1])
%     xlabel('Vs mm/sec')
%     ylabel('Vm (mV)')
%     
%     subplot(4,2,4);
%     plot(saveBinsFRVs{1}(:,1),saveBinsFRVs{1}(:,2), 'color', [0.5 0.5 0.5])
%     hold on
%     plot(saveBinsFRVs{2}(:,1),saveBinsFRVs{2}(:,2),'color', [0.5 0 0])
%     plot(saveBinsFRVs{3}(:,1),saveBinsFRVs{3}(:,2),'color', [1 0 0])
% 
%     %xlim([minValPlot+1 maxValPlot-1])
%     xlabel('Vs mm/sec')
%     ylabel('spikes/sec')
    
% yaw velocity 
figure();clf;
    subplot(2,1,1);
    plot(saveBinsVmVy{1}(:,1),saveBinsVmVy{1}(:,2), 'color', [0.5 0.5 0.5])  
    hold on 
    plot(saveBinsVmVy{2}(:,1),saveBinsVmVy{2}(:,2),'color', [0.5 0 0])
    plot(saveBinsVmVy{3}(:,1),saveBinsVmVy{3}(:,2),'color', [1 0 0])

    %xlim([minValPlot+1 maxValPlot-1])
    xlabel('Vy mm/sec')
    ylabel('Vm (mV)')
    set(gcf,'color',[1 1 1])
    box off
    
    subplot(2,1,2);
    plot(saveBinsFRVy{1}(:,1),saveBinsFRVy{1}(:,2), 'color', [0.5 0.5 0.5])
    hold on
    plot(saveBinsFRVy{2}(:,1),saveBinsFRVy{2}(:,2),'color', [0.5 0 0])
    plot(saveBinsFRVy{3}(:,1),saveBinsFRVy{3}(:,2),'color', [1 0 0])

    %xlim([minValPlot+1 maxValPlot-1])
    xlabel('Vy mm/sec')
    ylabel('spikes/sec')
    set(gcf,'color',[1 1 1])
    box off
    
    if savePlot
        saveas(gcf, fullfile(folder,'figures','diffHeadsActivityvsBehaviour_lineplots.fig'));
        saveas(gcf, fullfile(folder,'figures','diffHeadsActivityvsBehaviour_lineplots.svg'));
    end
    
% %% plot activity velocity relationships without splitting up b/w headings
% 
% 
    saveBinsVmSpeed_nohead = [];
    saveBinsFRSpeed_nohead = [];
    saveBinsVmVf_nohead = [];
    saveBinsFRVf_nohead = [];
    saveBinsVmVs_nohead = [];
    saveBinsFRVs_nohead = [];
    saveBinsVmVy_nohead = [];
    saveBinsFRVy_nohead = [];

edges = [-0.5:1:10];
    [N, edges, bin] = histcounts(speed, edges);
    tempFR = accumarray(bin+1, fRate_sec, [length(edges) 1]);
    tempVm = accumarray(bin+1, Vm, [length(edges) 1]);
    mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
    mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
    centers = edges(1:end-1)+diff(edges)/2;

    saveBinsVmSpeed_nohead = [centers' mean_binVm];
    saveBinsFRSpeed_nohead = [centers' mean_binFR];

edges = [0:0.5:8]; 
    [N, edges, bin] = histcounts(Vs(:,1), edges);
    tempFR = accumarray(bin+1, fRate_sec, [length(edges) 1]);
    tempVm = accumarray(bin+1, Vm, [length(edges) 1]);
    mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
    mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
    centers = edges(1:end-1)+diff(edges)/2;

    saveBinsVmVs_nohead = [centers' mean_binVm];
    saveBinsFRVs_nohead = [centers' mean_binFR];
                
edges = [-4:1:10]; 
    [N, edges, bin] = histcounts(Vf, edges);
    tempFR = accumarray(bin+1, fRate_sec, [length(edges) 1]);
    tempVm = accumarray(bin+1, Vm, [length(edges) 1]);
    mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
    mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
    centers = edges(1:end-1)+diff(edges)/2;

    saveBinsVmVf_nohead = [centers' mean_binVm];
    saveBinsFRVf_nohead = [centers' mean_binFR];
    
edges = [0:5:200]; 
    [N, edges, bin] = histcounts(Vy(:,1), edges);
    tempFR = accumarray(bin+1, fRate_sec, [length(edges) 1]);
    tempVm = accumarray(bin+1, Vm, [length(edges) 1]);
    mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
    mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
    centers = edges(1:end-1)+diff(edges)/2;

    saveBinsVmVy_nohead = [centers' mean_binVm];
    saveBinsFRVy_nohead = [centers' mean_binFR];

%%    
k = figure(); clf; 
    subplot(4,2,1);
    plot(saveBinsVmSpeed_nohead(:,1),saveBinsVmSpeed_nohead(:,2))
    xlabel('Speed mm/sec')
    ylabel('Vm (mV)')
    xlim([0 inf])
    %title(strcat(name, ' all headings'));

    subplot(4,2,3);
    plot(saveBinsFRSpeed_nohead(:,1),saveBinsFRSpeed_nohead(:,2))
    xlabel('Speed mm/sec')
    ylabel('spikes/sec')
    xlim([0 inf])
    
    subplot(4,2,5);
    plot(saveBinsVmVf_nohead(:,1),saveBinsVmVf_nohead(:,2))
    xlabel('Vf mm/sec')
    ylabel('Vm (mV)')

    subplot(4,2,7);
    plot(saveBinsFRVf_nohead(:,1),saveBinsFRVf_nohead(:,2))
    xlabel('Vf mm/sec')
    ylabel('spikes/sec')
    
    subplot(4,2,2);
    plot(saveBinsVmVs_nohead(:,1),saveBinsVmVs_nohead(:,2))
    xlabel('Vs mm/sec')
    ylabel('Vm (mV)')
    
    subplot(4,2,4);
    plot(saveBinsFRVs_nohead(:,1),saveBinsFRVs_nohead(:,2))
    xlabel('Vs mm/sec')
    ylabel('spikes/sec')
    
    subplot(4,2,6);
    plot(saveBinsVmVy_nohead(:,1),saveBinsVmVy_nohead(:,2))
    xlabel('Vy deg/sec')
    ylabel('Vm (mV)')

    subplot(4,2,8);
    plot(saveBinsFRVy_nohead(:,1),saveBinsFRVy_nohead(:,2))
    xlabel('Vy deg/sec')
    ylabel('spikes/sec')
    
    if savePlot
        saveas(gcf, fullfile(folder,'figures','activityvsBehaviour_all.fig'));
        saveas(gcf, fullfile(folder,'figures','activityvsBehaviour_all.svg'));
    end

end
