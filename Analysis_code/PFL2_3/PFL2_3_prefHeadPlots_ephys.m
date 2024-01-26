function h = PFL2_3_prefHeadPlots_ephys(rootDir, savePlots)
folders = get_folders_ephys(rootDir);

for f = 1:size(folders,1)

    folder = folders(f).folder; 
    if strcmp(folder(end),'.')
        folder = folder(1:end-2); 
    end
    
    processedDir = fullfile(folder,'processedData');

    load(fullfile(processedDir,'pro_behaviourData.mat'),'pro_behaviourData')
    load(fullfile(processedDir,'pro_trialData.mat'),'pro_trialData')
    
    bData = pro_behaviourData{1};
    tData = pro_trialData{1};

    edges_angle = [-180:30:180]; 
try
    total_mov_mm = abs(bData.vel_for) + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw))*4.5;
    no0vel_idx = find(total_mov_mm > 2); 
    angle = bData.angle; 
    angle = angle(no0vel_idx); 
        
    sum_mean_fRate{1} = zeros(length(edges_angle)-1,1);

    activity = tData.fRate_sec;
    activity = activity(no0vel_idx);   


    %angle
    [zscore, centers_angle, ~] = binData(activity, angle, edges_angle);
    sum_mean_fRate{1} = zscore; 
    smooth_FR = smoothdata(zscore,'gaussian',3);

    [bump_FR, model_data_FR, x_grid_FR] = fit_sinusoid(smooth_FR, [0,2*pi], 0);
    [bump_FR_vM, model_data_FR_vM, x_grid_FR_vM] = fit_von_Mises(smooth_FR, [0,2*pi], 0);
    
%     figure();
%     plot(centers_angle,zscore)
%     hold on
%     plot(centers_angle,smooth_FR)
%     plot(rad2deg((x_grid_FR))-165, feval(model_data_FR, x_grid_FR)); 
%     plot(rad2deg((x_grid_FR_vM))-165, feval(model_data_FR_vM, x_grid_FR_vM)); 
    
    sum_mean_Vm{1} = zeros(length(edges_angle)-1,1);
    
    try
        activity = tData.smooth_Vm;
    catch
        activity = tData.smoothVm;
    end
    activity = activity(no0vel_idx);   


    %angle
    [zscore, centers_angle, ~] = binData(activity, angle, edges_angle);
    sum_mean_Vm{1} = zscore;
    smooth_Vm = smoothdata(zscore,'gaussian',3);
    [bump_Vm, model_data_Vm, x_grid_Vm] = fit_sinusoid(smooth_Vm, [0,2*pi], 0);
    [bump_Vm_vM, model_data_Vm_vM, x_grid_Vm_vM] = fit_von_Mises(smooth_Vm, [0,2*pi], 0);
%     figure();
%     plot(centers_angle,zscore)
%     hold on
%     plot(centers_angle,smooth_Vm)
%     plot(rad2deg((x_grid_Vm))-165, feval(model_data_Vm, x_grid_Vm));
%     plot(rad2deg((x_grid_Vm_vM))-165, feval(model_data_Vm_vM, x_grid_Vm_vM)); 
    


    figure();
    set(gcf,'color','w')
    set(gcf,'Renderer','painters')
    
%     subplot(2,3,1)
%     plot(centers_angle,smooth_FR)
%     xlabel('cue pos')
%     ylabel('FR')
%     box off
%     subplot(2,3,2)
%     plot(rad2deg((x_grid_FR))-165, feval(model_data_FR, x_grid_FR));  
%     xlabel('cue pos')
%     ylabel('FR_vonMisesFit')
%     annotation('textbox', [0.55, 0.8, 0.1, 0.1], 'String', num2str(bump_FR.adj_rs),'EdgeColor','none')
%     box off
%     subplot(2,3,3)
%     plot(rad2deg((x_grid_FR_vM))-165, feval(model_data_FR_vM, x_grid_FR_vM));
%     xlabel('cue pos')
%     ylabel('FR_sinusoidFit')
%     annotation('textbox', [0.85, 0.8, 0.1, 0.1], 'String', num2str(bump_FR_vM.adj_rs),'EdgeColor','none')
%     box off
    subplot(2,1,1)
    plot(centers_angle,smooth_Vm)
    xlabel('cue pos')
    ylabel('Vm')
    box off
    subplot(2,1,2)
    plot(rad2deg((x_grid_Vm))-165, feval(model_data_Vm, x_grid_Vm));  
    xlabel('cue pos')
    ylabel('sinusoidFit')
    annotation('textbox', [0.8, 0.3, 0.1, 0.1], 'String', num2str(bump_Vm.rs),'EdgeColor','none')
    box off
%     subplot(2,3,6)
%     plot(rad2deg((x_grid_Vm_vM))-165, feval(model_data_Vm_vM, x_grid_Vm_vM)); 
%     xlabel('cue pos')
%     ylabel('Vm_sinusoidFit')
%     annotation('textbox', [0.88, 0.4, 0.1, 0.1], 'String', num2str(bump_Vm_vM.adj_rs),'EdgeColor','none')
%     box off


    if savePlots
        saveas(gcf,fullfile(folder,'figures','prefHeadPlots.fig'))
    end
catch
    disp([folder,' failed'])
end
end
end