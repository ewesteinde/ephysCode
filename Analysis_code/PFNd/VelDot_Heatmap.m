function VelDot_Heatmap(tStart, tEnd, processed_behaviourData, processed_trialData, fileName)


for type = 1:2
    activity_values = [];
    x_values = [];
    y_values = [];
    for t = tStart:tEnd
        
        tData = processed_trialData{t};
    
        if type == 1
            data = getfield(tData,'fRate_sec');
        else
            data = getfield(tData,'smooth_Vm');
        end
        
        activity_values = cat(1,activity_values, data);
        x_values = cat(1, x_values, (processed_behaviourData{t}.angle));
        y_values = cat(1, y_values, processed_behaviourData{t}.vel_dot);

    end

    x_edges = [-180:10:180]; %
    y_edges = [-2:.5:10];


    if type == 1
        [N_fRate, heatmap_fRate, x_centers_fRate, y_centers_fRate] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    else
        [N_Vm, heatmap_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    end
 
end


 figure(); clf;
            f(1) = subplot(2,1,1);
            imagesc(flip(heatmap_fRate))
            ylabel('Vdot mm/s')
            xt = get(gca, 'XTick');
            xtnew = linspace(min(xt), max(xt), 7);                             
            xtlbl = linspace(-135, 165, numel(xtnew));                  
            set(gca, 'XTick',xtnew, 'XTickLabel',xtlbl)
            xlabel('angle')
            colorbar
            title('fRate')
            yt = get(gca, 'YTick');
            ytnew = linspace(min(yt),max(yt),6);
            ytlbl = linspace(8.5, -1.5, numel(ytnew));
            set(gca, 'YTick',ytnew, 'YTickLabel',ytlbl)
            
            f(2) = subplot(2,1,2);
            imagesc(flip(heatmap_Vm))
            ylabel('Vdot mm/s')
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
            
keep = input('Save? ','s');

if strcmp(keep, 'y')
    cd(fileName) 
    saveas(gcf,'directionPref_velDot.fig')
    cd('/Users/elenawesteinde/Documents/EphysCode') %'/Users/elenawesteinde/Documents/EphysCode/Analysis_code' 'C:\Code\EphysCode'
end
end