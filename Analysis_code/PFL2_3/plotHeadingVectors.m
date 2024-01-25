function plotHeadingVectors(folder,processed_behaviourData, trial, sampRate, minVel,savePlot)

%% overlapping 60 second windows slid by ~1s increments
    window = 60 * sampRate; 
    mean_headingVectors = [];
    count = 1; 
    speed = sqrt(processed_behaviourData.vel_for.^2 + processed_behaviourData.vel_side.^2);
    for i = 1:sampRate:length(processed_behaviourData.angle)
        idx = round(i - window/2:i + window/2); 
        if idx(end) > length(processed_behaviourData.angle)
            idx = idx(1):1:length(processed_behaviourData.angle);
        end
        if idx(1) < 1 
            idx = 1:idx(end); 
        end
        angle_temp = processed_behaviourData.angle(idx); 
        speed_temp = speed(idx); 
        angles_flyFor = angle_temp(speed_temp > minVel); 
        if ~isempty(angles_flyFor)
            x = cosd(angles_flyFor); 
            y = sind(angles_flyFor); %my arena has - angles to the left of the fly, + to the right, multiply y component by -1 to align physical arena coords to polar plot angles
            mean_headingVectors(1,count)= sum(x)/length(x); 
            mean_headingVectors(2,count)= sum(y)/length(y); 
            count = count + 1; 
        end
    end 

            trial_angle = processed_behaviourData.angle(speed > minVel);
            xTrial = cosd(trial_angle);
            yTrial = sind(trial_angle); 
            trialHeadingVector(1) = sum(xTrial)/length(xTrial);
            trialHeadingVector(2) = sum(yTrial)/length(yTrial);

            rhoTrial = sqrt(trialHeadingVector(1).^2 + trialHeadingVector(2).^2); 
            thetaTrial = atan2(trialHeadingVector(2),trialHeadingVector(1)); 


            rho = sqrt(mean_headingVectors(1,:).^2 + mean_headingVectors(2,:).^2); 
            theta = atan2(mean_headingVectors(2,:),mean_headingVectors(1,:)); 

%             rho_all = [rho_all rho];
%             theta_all = [theta_all theta];

% uncomment to plot each trial's heading data one by one
%             polarscatter(theta, rho,'o')
%             %rlim([0 1])
%             hold on
%             colormap(gca, 'hot')

                figure();clf; 
                title([folder(end-10:end), 'trial ',num2str(trial)], 'Interpreter', 'none')
                polarscatter(theta, rho,'o')
                rlim([0 1])
                hold on
                polarplot([0,thetaTrial], [0,rhoTrial],'k','LineWidth',2)
                str = {['mean vector str: ',num2str(rhoTrial)]};
                dim = [0.7 0.72 0.2 0.2];
                annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','w')
                set(gcf,'color','w');
                set(gcf,'renderer','painters')
                rlim([0 1])
                
                
                if savePlot
                    saveas(gcf, fullfile(folder,'figures',['trialVectors_trial_',num2str(trial),'.fig']));
                    saveas(gcf, fullfile(folder,'figures',['trialVectors_trial_',num2str(trial),'.svg']));
                end
end

