function  [thetaTrial, rhoTrial, theta, rho] = plotHeading(behaviourData, windowSize, sampRate, folder, plot) 

    if iscell(behaviourData)
        behaviourData = behaviourData{1};
    end

    window = windowSize * round(sampRate); 
    if window > length(behaviourData.vel_for)
        window = length(behaviourData.vel_for); 
    end
    count = 1; 

    %% overlapping 60 second windows slid by ~1s increments
    mean_headingVectors = [];
    idx_windows = [];
    speed = sqrt(behaviourData.vel_for.^2 + behaviourData.vel_side.^2);
    for i = 1:round(sampRate):length(behaviourData.angle) - window + 1
        idx = i - window/2:i + window/2; 
        if idx(end) > length(behaviourData.angle)
            idx = idx(1):1:length(behaviourData.angle);
        elseif idx(1) < 1 
            idx = 1:idx(end); 
        end
        angle_temp = behaviourData.angle(idx); 
        speed_temp = speed(idx); 
        angles_flyFor = angle_temp(speed_temp > 1.5); 
        if ~isempty(angles_flyFor)
            x = cosd(angles_flyFor); 
            y = sind(angles_flyFor); %my arena has - angles to the left of the fly, + to the right, multiply y component by -1 to align physical arena coords to polar plot angles
            idx_windows(count,1) = idx(1);
            idx_windows(count,2) = idx(end);
            idx_windows(count,3) = length(find(speed_temp > 1.5))/length(speed_temp)*100;
            mean_headingVectors(1,count)= sum(x)/length(x); 
            mean_headingVectors(2,count)= sum(y)/length(y); 
            count = count + 1; 
        end
    end 

    trial_angle = behaviourData.angle(speed > 1.5);
    xTrial = cosd(trial_angle);
    yTrial = sind(trial_angle); 
    trialHeadingVector(1) = sum(xTrial)/length(xTrial);
    trialHeadingVector(2) = sum(yTrial)/length(yTrial);

    rhoTrial = sqrt(trialHeadingVector(1).^2 + trialHeadingVector(2).^2); 
    thetaTrial = atan2(trialHeadingVector(2),trialHeadingVector(1)); 


    rho = sqrt(mean_headingVectors(1,:).^2 + mean_headingVectors(2,:).^2); 
    theta = atan2(mean_headingVectors(2,:),mean_headingVectors(1,:)); 

    figure();clf; 
    polarscatter(theta, rho,'o')
    rlim([0 1])
    %colormap(gca, 'hot')
    hold on
    polarplot([0,thetaTrial], [0,rhoTrial],'k','LineWidth',2)
    str = {['mean vector str: ',num2str(rhoTrial)]};
    dim = [0.7 0.72 0.2 0.2];
    annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','w')
    set(gcf,'color','w');
    rlim([0 1])

    if plot
        saveas(gcf, fullfile(folder,'figures','headingDist_trialVector.fig'));
        saveas(gcf, fullfile(folder,'figures','headingDist_trialVector.svg'));
    end
end