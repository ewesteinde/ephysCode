function WholeTrialFig(processed_behaviourData,processed_trialData,figure_dir,trial,savePlot)

    figure();clf;  
    set(gcf,'renderer','painters')
    set(gcf,'color','w')
 
    h(1) = subplot(3,1,1);
    plot(processed_behaviourData.time, processed_behaviourData.angle, 'k') 
    ylabel('angle')
    ylim([-180 180])
    box off
    title([figure_dir(end-17:end-7), 'trial ',num2str(trial)], 'Interpreter', 'none')
    
    h(2) = subplot(3,1,2);
    yyaxis right
    plot(processed_behaviourData.time, processed_behaviourData.vel_for) 
    ylabel('for (mm/s)')
    yyaxis left 
    plot(processed_behaviourData.time, processed_behaviourData.vel_yaw) 
    ylabel('yaw (deg/s)')
    %ylim([-8 8])
    box off


    h(3) = subplot(3,1,3);
    plot(processed_behaviourData.time, processed_trialData.scaledOutput, 'k') 
    ylabel('Vm (mV)')
    set(gcf,'color','w')
    box off

    linkaxes(h,'x');
    
    if savePlot
        saveas(gcf, fullfile(figure_dir,['whole_trial_',num2str(trial),'.fig']))
    end
end