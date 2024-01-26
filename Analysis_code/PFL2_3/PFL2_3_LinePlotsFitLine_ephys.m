function [vfcoeff,vycoeff] = PFL2_3_LinePlotsFitLine_ephys(folder, bData, tData, savePlots)

    edges_vy = [-200:20:200];
    edges_vf = [-4:0.5:10];

    total_mov_mm = abs(bData.vel_for) + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw))*4.5;
    no0vel_idx = find(total_mov_mm > 2);
    vf = bData.vel_for;
    vf = vf(no0vel_idx); 
    vy = bData.vel_yaw;
    vy = vy(no0vel_idx); 
  
    sum_mean_Vm{1} = zeros(length(edges_vy)-1,1);
    sum_mean_Vm{2} = zeros(length(edges_vf)-1,1);

    try
        activity = tData.smooth_Vm;
    catch
        activity = tData.smoothVm;
    end
    activity = activity(no0vel_idx);   


    %vf 
    behaviour = vf; 
    [zscore, centers_vf, ~] = binData(activity, behaviour, edges_vf);
    sum_mean_Vm{2} = zscore;

    %vy 
     behaviour = vy; 
    [zscore, centers_vy, ~] = binData(activity, behaviour, edges_vy);
    sum_mean_Vm{1} = zscore; 
    

    figure();
    set(gcf,'color','w')
    set(gcf,'Renderer','painters')
    subplot(2,1,1)
    hold on
    plot(centers_vf,sum_mean_Vm{2},'k')
    vfcoeff = polyfit(centers_vf(~isnan(sum_mean_Vm{2})), sum_mean_Vm{2}(~isnan(sum_mean_Vm{2})), 1);
    xFit = linspace(min(centers_vf(~isnan(sum_mean_Vm{2}))), max(centers_vf(~isnan(sum_mean_Vm{2}))), 1000);
    yFit = polyval(vfcoeff, xFit);
    plot(xFit,yFit,'k')
    xlabel('vf mm/s')
    ylabel('Vm')
    box off
    
    subplot(2,1,2)
    hold on
    plot(centers_vy,sum_mean_Vm{1},'k')
    vycoeff = polyfit(centers_vy(~isnan(sum_mean_Vm{1})), sum_mean_Vm{1}(~isnan(sum_mean_Vm{1})), 1);
    xFit = linspace(min(centers_vy(~isnan(sum_mean_Vm{1}))), max(centers_vy(~isnan(sum_mean_Vm{1}))), 1000);
    yFit = polyval(vycoeff , xFit);
    plot(xFit,yFit,'k')
    xlabel('vy deg/s')
    ylabel('Vm')
    box off
    


    if savePlots
        saveas(gcf,fullfile(folder,'figures','LinePlotsLineFit.fig'))
    end
end