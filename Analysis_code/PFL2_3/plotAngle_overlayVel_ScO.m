function plotAngle_overlayVel_ScO(processed_trialData, processed_behaviourData)

nActivity = processed_trialData.scaledOutput_down;
figure();clf;      
        
        h(1) = subplot(4,1,1);
        yyaxis left 
        plot(processed_behaviourData.time,nActivity, 'k') 
        yyaxis right
        plot(processed_behaviourData.time, processed_behaviourData.angle, 'r') 
        ylabel('angle')
       
        
        h(2) = subplot(4,1,2);
        yyaxis left
        plot(processed_behaviourData.time,nActivity, 'k') 
        ylabel('Vm')
        %ylim([-70 -45])
        hold on 
        yyaxis right
        plot(processed_behaviourData.time, processed_behaviourData.vel_for, 'r')
        ylim([-(max(processed_behaviourData.vel_for)) max(processed_behaviourData.vel_for)])
        ylabel('Vf mm/sec')
        
        
        
        h(3) = subplot(4,1,3);
        yyaxis left
        plot(processed_behaviourData.time,nActivity, 'k')
        ylabel('Vm')
       % ylim([-70 -45])
        hold on 
        yyaxis right
        plot(processed_behaviourData.time, processed_behaviourData.vel_side, 'r')
        ylim([-(max(processed_behaviourData.vel_side)) max(processed_behaviourData.vel_side)])
        ylabel('Vs mm/sec')
        
        h(4) = subplot(4,1,4);
        yyaxis left
        plot(processed_behaviourData.time,nActivity, 'k')
        ylabel('Vm')

        %ylim([-55 -])
        hold on 
        yyaxis right
        plot(processed_behaviourData.time, processed_behaviourData.vel_yaw, 'r')
        ylim([(min(processed_behaviourData.vel_yaw)) max(processed_behaviourData.vel_yaw)])
        ylabel('Vy deg/sec')
        
        
        linkaxes(h,'x');
end