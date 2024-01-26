function [h,j] = Activity_Heatmap_ephys(folder, angleBin, processed_behaviourData, processed_trialData, savePlots)

tData = processed_trialData;
for type = 1:2

    if type == 1
        data = tData.fRate_sec;
    else
        try
            data = tData.smooth_Vm;
        catch
            data = tData.smoothVm;
        end
    end

    activity_values = data;
    x_values = processed_behaviourData.angle;
    y_values = processed_behaviourData.vel_for;

    

    x_edges = [-180:angleBin:180]; %
    y_edges = [-4:.5:10];


    if type == 1
        [N_fRate, heatmap_fRate, x_centers_fRate, y_centers_fRate] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    else
        [N_Vm, heatmap_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    end
 
end

heatmap_fRate(N_fRate == 0) = nan; 
heatmap_Vm(N_Vm == 0) = nan; 
ncol = length(unique(heatmap_fRate(~isnan(heatmap_fRate))));
color = flipud(cbrewer2('RdYlBu', ncol));

 
 h = figure();
            set(gcf,'color','w')
            set(gcf,'renderer','painters')
            subplot(2,2,1);
            s = pcolor(heatmap_fRate);
            colormap(color)
            hold on
            
            xt = linspace(1,numel(x_centers_fRate),5);                            
            xtlbl = linspace(x_centers_fRate(xt(1)), x_centers_fRate(xt(end)), 5);   
            set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
            ylabel('vf mm/s')
            xlabel('cue pos')
            colorbar
            yt = linspace(1,numel(y_centers_fRate),5); 
            ytlbl = linspace(y_centers_fRate(yt(1)), y_centers_fRate(yt(end)), 5);
            set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
            set(s, 'EdgeColor', 'none');
            set(gca,'color','none')
            box off
            
            subplot(2,2,3);
            ncol = length(unique(heatmap_Vm(~isnan(heatmap_Vm)))); 
            color = flipud(cbrewer2('RdYlBu', ncol));
            s = pcolor(heatmap_Vm);
            colormap(color)
            hold on
            
            xt = linspace(1,numel(x_centers_Vm),5);                            
            xtlbl = linspace(x_centers_Vm(xt(1)), x_centers_Vm(xt(end)), 5);   
            set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
            ylabel('vf mm/s')
            xlabel('cue pos')
            colorbar
            yt = linspace(1,numel(y_centers_Vm),5); 
            ytlbl = linspace(y_centers_Vm(yt(1)), y_centers_Vm(yt(end)), 5);
            set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
            set(s, 'EdgeColor', 'none');
            set(gca,'color','none')
            box off

for type = 1:2
    if type == 1
        data = tData.fRate_sec;
    else
        try
            data = tData.smooth_Vm;
        catch
            data = tData.smoothVm;
        end
    end

    activity_values = data;
    x_values = processed_behaviourData.vel_side;
    y_values = processed_behaviourData.vel_for;


    x_edges = [-10:.5:10]; 
    y_edges = [-4:0.5:10];


    if type == 1
        [N_fRate, heatmap_fRate, x_centers_fRate, y_centers_fRate] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    else
        [N_Vm, heatmap_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    end
 
end

heatmap_fRate(N_fRate == 0) = nan; 
heatmap_Vm(N_Vm == 0) = nan; 

heatmap_fRate(N_fRate == 0) = nan; 
heatmap_Vm(N_Vm == 0) = nan; 
ncol = length(unique(heatmap_fRate(~isnan(heatmap_fRate))));
color = flipud(cbrewer2('RdYlBu', ncol));

    subplot(2,2,2);
    s = pcolor(heatmap_fRate);
    colormap(color)
    hold on

    xt = linspace(1,numel(x_centers_fRate),7);                            
    xtlbl = linspace(x_centers_fRate(xt(1)), x_centers_fRate(xt(end)), 7);   
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
    ylabel('vf mm/s')
    xlabel('vs mm/s')
    colorbar
    yt = linspace(1,numel(y_centers_fRate),5); 
    ytlbl = linspace(y_centers_fRate(yt(1)), y_centers_fRate(yt(end)), 5);
    set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
    set(s, 'EdgeColor', 'none');
    set(gca,'color','none')
    box off

    subplot(2,2,4);
    ncol = length(unique(heatmap_Vm(~isnan(heatmap_Vm)))); 
    color = flipud(cbrewer2('RdYlBu', ncol));
    s = pcolor(heatmap_Vm);
    colormap(color)
    hold on

    xt = linspace(1,numel(x_centers_Vm),7);                            
    xtlbl = linspace(x_centers_Vm(xt(1)), x_centers_Vm(xt(end)), 7);   
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
    ylabel('vf mm/s')
    xlabel('vs mm/s')
    colorbar
    yt = linspace(1,numel(y_centers_Vm),5); 
    ytlbl = linspace(y_centers_Vm(yt(1)), y_centers_Vm(yt(end)), 5);
    set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
    set(s, 'EdgeColor', 'none');
    set(gca,'color','none')
    box off
    
    if savePlots
        saveas(gcf,fullfile(folder,'figures','heatMaps_vf_s.fig'))
    end

% uncomment to check yaw relationship


for type = 1:2
        if type == 1
        data = tData.fRate_sec;
        else
            try
                data = tData.smooth_Vm;
            catch
                data = tData.smoothVm;
            end
        end
        
        activity_values = data;
        x_values = processed_behaviourData.vel_yaw;
        y_values = processed_behaviourData.vel_for;


    x_edges = [-200:10:200]; 
    y_edges = [-4:0.5:10];


    if type == 1
        [N_fRate, heatmap_fRate, x_centers_fRate, y_centers_fRate] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    else
        [N_Vm, heatmap_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    end
 
end

j = figure();
            set(gcf,'color','w')
            set(gcf,'renderer','painters')


heatmap_fRate(N_fRate == 0) = nan; 
heatmap_Vm(N_Vm == 0) = nan;             
ncol = length(unique(heatmap_fRate(~isnan(heatmap_fRate)))); 
color = flipud(cbrewer2('RdYlBu', ncol));

    subplot(2,2,1);
    s = pcolor(heatmap_fRate);
    colormap(color)
    hold on

    xt = linspace(1,numel(x_centers_fRate),7);                            
    xtlbl = linspace(x_centers_fRate(xt(1)), x_centers_fRate(xt(end)), 7);   
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
    ylabel('vf mm/s')
    xlabel('vy deg/s')
    colorbar
    yt = linspace(1,numel(y_centers_fRate),5); 
    ytlbl = linspace(y_centers_fRate(yt(1)), y_centers_fRate(yt(end)), 5);
    set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
    set(s, 'EdgeColor', 'none');
    set(gca,'color','none')
    box off

    subplot(2,2,3);
    ncol = length(unique(heatmap_Vm(~isnan(heatmap_Vm)))); 
    color = flipud(cbrewer2('RdYlBu', ncol));
    s = pcolor(heatmap_Vm);
    colormap(color)
    hold on

    xt = linspace(1,numel(x_centers_Vm),7);                            
    xtlbl = linspace(x_centers_Vm(xt(1)), x_centers_Vm(xt(end)), 7);   
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
    ylabel('vf mm/s')
    xlabel('vy deg/s')
    colorbar
    yt = linspace(1,numel(y_centers_Vm),5); 
    ytlbl = linspace(y_centers_Vm(yt(1)), y_centers_Vm(yt(end)), 5);
    set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
    set(s, 'EdgeColor', 'none');
    set(gca,'color','none')
    box off

    

for type = 1:2

        if type == 1
        data = tData.fRate_sec;
        else
            try
                data = tData.smooth_Vm;
            catch
                data = tData.smoothVm;
            end
        end
        
        activity_values = data;
        x_values = processed_behaviourData.angle;
        y_values = processed_behaviourData.vel_yaw;

    x_edges = [-180:angleBin:180]; 
    y_edges = [-200:10:200];


    if type == 1
        [N_fRate, heatmap_fRate, x_centers_fRate, y_centers_fRate] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    else
        [N_Vm, heatmap_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    end
    
end

heatmap_fRate(N_fRate == 0) = nan; 
heatmap_Vm(N_Vm == 0) = nan;  
ncol = length(unique(heatmap_fRate(~isnan(heatmap_fRate))));
color = flipud(cbrewer2('RdYlBu', ncol));

    subplot(2,2,2);
    s = pcolor(heatmap_fRate);
    colormap(color)
    hold on

    xt = linspace(1,numel(x_centers_fRate),5);                            
    xtlbl = linspace(x_centers_fRate(xt(1)), x_centers_fRate(xt(end)), 5);   
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
    xlabel('cue pos')
    ylabel('vy deg/s')
    colorbar
    yt = linspace(1,numel(y_centers_fRate),5); 
    ytlbl = linspace(y_centers_fRate(yt(1)), y_centers_fRate(yt(end)), 5);
    set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
    set(s, 'EdgeColor', 'none');
    set(gca,'color','none')
    box off

    subplot(2,2,4);
    ncol = length(unique(heatmap_Vm(~isnan(heatmap_Vm)))); 
    color = flipud(cbrewer2('RdYlBu', ncol));
    s = pcolor(heatmap_Vm);
    colormap(color)
    hold on

    xt = linspace(1,numel(x_centers_Vm),5);                            
    xtlbl = linspace(x_centers_Vm(xt(1)), x_centers_Vm(xt(end)), 5);   
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
    xlabel('cue pos')
    ylabel('vy deg/s')
    colorbar
    yt = linspace(1,numel(y_centers_Vm),5); 
    ytlbl = linspace(y_centers_Vm(yt(1)), y_centers_Vm(yt(end)), 5);
    set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
    set(s, 'EdgeColor', 'none');
    set(gca,'color','none')
    box off

    if savePlots
        saveas(gcf,fullfile(folder,'figures','heatMaps_vy.fig'))
    end
end