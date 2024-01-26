function h = PFL2_3_basicLinePlots_ephys(folder, bData, tData, savePlots)

    edges_vy = [-200:10:200];
    edges_vs = [-8:0.5:8];
    edges_vf = [-4:0.5:10];
    edges_angle = [-180:20:180]; 

    total_mov_mm = abs(bData.vel_for) + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw))*4.5;
    no0vel_idx = find(total_mov_mm > 0);
    vf = bData.vel_for;
    vf = vf(no0vel_idx); 
    vs = bData.vel_side;
    vs = vs(no0vel_idx); 
    vy = bData.vel_yaw;
    vy = vy(no0vel_idx); 
    angle = bData.angle; 
    angle = angle(no0vel_idx); 
        
    sum_mean_fRate{1} = zeros(length(edges_vy)-1,1);
    sum_mean_fRate{2} = zeros(length(edges_angle)-1,1);
    sum_mean_fRate{3} = zeros(length(edges_vf)-1,1);
    sum_mean_fRate{4} = zeros(length(edges_vs)-1,1);

    activity = tData.fRate_sec;
    activity = activity(no0vel_idx);   


    %vf 
    behaviour = vf; 
    [zscore, ~, ~] = binData(activity, behaviour, edges_vf);
    sum_mean_fRate{3} = zscore;

    %vy 
     behaviour = vy; 
    [zscore, ~, ~] = binData(activity, behaviour, edges_vy);
    sum_mean_fRate{1} = zscore; 
    
    %vs
     behaviour = vs; 
    [zscore, ~, ~] = binData(activity, behaviour, edges_vs);
    sum_mean_fRate{4} = zscore; 

    %angle
    [zscore, ~, ~] = binData(activity, angle, edges_angle);
    sum_mean_fRate{2} = zscore; 
    
    sum_mean_Vm{1} = zeros(length(edges_vy)-1,1);
    sum_mean_Vm{2} = zeros(length(edges_angle)-1,1);
    sum_mean_Vm{3} = zeros(length(edges_vf)-1,1);
    sum_mean_Vm{4} = zeros(length(edges_vs)-1,1);

    try
        activity = tData.smooth_Vm;
    catch
        activity = tData.smoothVm;
    end
    activity = activity(no0vel_idx);   


    %vf 
    behaviour = vf; 
    [zscore, centers_vf, ~] = binData(activity, behaviour, edges_vf);
    sum_mean_Vm{3} = zscore;

    %vy 
     behaviour = vy; 
    [zscore, centers_vy, ~] = binData(activity, behaviour, edges_vy);
    sum_mean_Vm{1} = zscore; 
    
    %vs
     behaviour = vs; 
    [zscore, centers_vs, ~] = binData(activity, behaviour, edges_vs);
    sum_mean_Vm{4} = zscore; 

    %angle
    [zscore, centers_angle, ~] = binData(activity, angle, edges_angle);
    sum_mean_Vm{2} = zscore;

    h = figure();
    set(gcf,'color','w')
    set(gcf,'Renderer','painters')
    subplot(4,2,1)
    plot(centers_vf,sum_mean_fRate{3})
    xlabel('vf mm/s')
    ylabel('FR')
    box off
    subplot(4,2,2)
    plot(centers_vf,sum_mean_Vm{3})
    xlabel('vf mm/s')
    ylabel('Vm')
    box off
    
    subplot(4,2,3)
    plot(centers_vy,sum_mean_fRate{1})
    xlabel('vy deg/s')
    ylabel('FR')
    box off
    subplot(4,2,4)
    plot(centers_vy,sum_mean_Vm{1})
    xlabel('vy deg/s')
    ylabel('Vm')
    box off
    
    subplot(4,2,5)
    plot(centers_vs,sum_mean_fRate{4})
    xlabel('vs mm/s')
    ylabel('FR')
    box off
    subplot(4,2,6)
    plot(centers_vs,sum_mean_Vm{4})
    xlabel('vs mm/s')
    ylabel('Vm')
    box off
    
    subplot(4,2,7)
    plot(centers_angle,sum_mean_fRate{2})
    xlabel('cue pos')
    ylabel('FR')
    box off
    subplot(4,2,8)
    plot(centers_angle,sum_mean_Vm{2})
    xlabel('cue pos')
    ylabel('Vm')
    box off


    if savePlots
        saveas(gcf,fullfile(folder,'figures','basicLinePlots.fig'))
    end
end