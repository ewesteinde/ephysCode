  startOL = 178.189*1000; % sec
    endOL = 203.249*1000;
    lengthCL = 179.99*1000;
    
% start = 1819.45;
% finish = 1844.47;
%     startOL  = find(processed_behaviourData.time == start); 
%     endOL = find(processed_behaviourData.time == finish); 

[CL, OL, CL_startStopIdx, OL_startStopIdx] = separateCL_OL(startOL, endOL, lengthCL, processed_behaviourData, processed_trialData);
%%

for chunk = 1:length(OL_startStopIdx)

    start = OL_startStopIdx(chunk, 1);
    finish = OL_startStopIdx(chunk, 2);
    
    f = fieldnames(processed_behaviourData);
    for k=1:numel(f)
        chunk_behaviourData.(f{k}) = processed_behaviourData.(f{k})(start:finish); 
    end


    f = fieldnames(processed_trialData);
    for k=1:numel(f)
        chunk_trialData.(f{k}) = processed_trialData.(f{k})(start:finish); 
    end
    
    figure(3);
    subplot(length(OL_startStopIdx),1,chunk)
        Vm = chunk_trialData.smooth_Vm;
        angle = chunk_behaviourData.angle;
        edges = [-180:30:180];
        [centers, mean_bin] = create_binned_mean(angle, Vm, edges);
        bump_amp = max(mean_bin) - min(mean_bin);
        plot(centers, mean_bin,'-o');
        text(-175, -50, num2str(bump_amp))
        ylim([-64 -48])
        xlabel('Angle');
        ylabel('Vm'); 
        title('OL heading pref')
    figure(4);
        subplot(length(OL_startStopIdx),1,chunk)
  
        activity_values = chunk_trialData.smooth_Vm';
        x_values = chunk_behaviourData.angle';
        y_values = chunk_behaviourData.vel_for';
        x_edges = [-180:30:180]; %
        y_edges = [-2:.5:10];
        [N_Vm, heatmap_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
        imagesc(flip(heatmap_Vm))
            
        ylabel('Vf mm/s')
        xt = get(gca, 'XTick');
        xtnew = linspace(min(xt), max(xt), 7);                             
        xtlbl = linspace(-135, 165, numel(xtnew));                  
        set(gca, 'XTick',xtnew, 'XTickLabel',xtlbl)
        xlabel('angle')
        colorbar
        title('Vm')
        yt = get(gca, 'YTick');
        ytnew = linspace(min(yt),max(yt),6);
        ytlbl = linspace(8.5, -1.5, numel(ytnew));
        set(gca, 'YTick',ytnew, 'YTickLabel',ytlbl)
        set(gcf,'color','w')

end

for chunk = 1:length(CL_startStopIdx)

    start = CL_startStopIdx(chunk, 1);
    finish = CL_startStopIdx(chunk, 2);
    
    f = fieldnames(processed_behaviourData);
    for k=1:numel(f)
        chunk_behaviourData.(f{k}) = processed_behaviourData.(f{k})(start:finish); 
    end


    f = fieldnames(processed_trialData);
    for k=1:numel(f)
        chunk_trialData.(f{k}) = processed_trialData.(f{k})(start:finish); 
    end
    
    figure(5);
    subplot(length(CL_startStopIdx),1,chunk)
        Vm = chunk_trialData.smooth_Vm;
        angle = chunk_behaviourData.angle;
        edges = [-180:30:180];
        [centers, mean_bin] = create_binned_mean(angle, Vm, edges);
        bump_amp = max(mean_bin) - min(mean_bin);
        plot(centers, mean_bin,'-o');
        text(-175, -50, num2str(bump_amp))
        ylim([-64 -48])
        xlabel('Angle');
        ylabel('Vm'); 
        title('CL heading pref')
    figure(6);
        subplot(length(CL_startStopIdx),1,chunk)
  
        activity_values = chunk_trialData.smooth_Vm';
        x_values = chunk_behaviourData.angle';
        y_values = chunk_behaviourData.vel_for';
        x_edges = [-180:30:180]; %
        y_edges = [-2:.5:10];
        [N_Vm, heatmap_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
        imagesc(flip(heatmap_Vm))
            
        ylabel('Vf mm/s')
        xt = get(gca, 'XTick');
        xtnew = linspace(min(xt), max(xt), 7);                             
        xtlbl = linspace(-135, 165, numel(xtnew));                  
        set(gca, 'XTick',xtnew, 'XTickLabel',xtlbl)
        xlabel('angle')
        colorbar
        title('Vm')
        yt = get(gca, 'YTick');
        ytnew = linspace(min(yt),max(yt),6);
        ytlbl = linspace(8.5, -1.5, numel(ytnew));
        set(gca, 'YTick',ytnew, 'YTickLabel',ytlbl)
        set(gcf,'color','w')

end
   