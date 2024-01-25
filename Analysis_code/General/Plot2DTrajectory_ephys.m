function [xPos,yPos] = Plot2DTrajectory_ephys(behaviourData,folder,s)
    newSampRate = 60; % check
    
    vq = behaviourData.time(1):1/newSampRate:behaviourData.time(end);
    fn = fieldnames(behaviourData);
    if istable(behaviourData)
        fn = fn(1:end - 3); 
    end
    for field = 1:numel(fn)
        if isempty(regexp(fn{field},'disp'))
            bData_interp.(fn{field}) = interp1(behaviourData.time,behaviourData.(fn{field}),vq,'spline');
        end
    end  
    

    yawAngPos = bData_interp.angle;
    fwdAngVel = bData_interp.vel_for;
    slideAngVel = bData_interp.vel_side;
    if s
        stim = bData_interp.stim;
    end


    % conversion factor between degrees and mm
    circum = 9 * pi; % circumference of ball, in mm
    mmPerDeg = circum / 360; % mm per degree of ball

    % position incorporating heading - as if fly were walking on x-y plane,
    %  x-y coordinates at each time point
    % start with fly at (0,0) and facing 0 deg
    zeroedYawAngPos = yawAngPos - yawAngPos(1); 

    % movement in x (in degrees) at each time point
    xChangePos = (fwdAngVel ./ newSampRate) .* sind(zeroedYawAngPos) + ...
        (slideAngVel ./ newSampRate) .* sind(zeroedYawAngPos + 90);  

    % x position in mm (i.e. x-coordinate of fly's position at each time 
    %  point), starts at 0
    xPos = (cumsum(xChangePos) - xChangePos(1)) .* mmPerDeg;
    minX = min(xPos);
    maxX = max(xPos);

    % movement in y (in degrees) at each time point
    yChangePos = (fwdAngVel ./ newSampRate) .* cosd(zeroedYawAngPos) + ...
        (slideAngVel ./ newSampRate) .* cosd(zeroedYawAngPos + 90);

    % y position in mm (i.e. y-coordinate of fly's position at each time 
    %  point), starts at 0
    yPos = (cumsum(yChangePos) - yChangePos(1)) .* mmPerDeg;
    minY = min(yPos);
    maxY = max(yPos);

    %[~, ~, transition] = detect_jumps(ftT, 2, 2);
    %time = seconds(ftT.trialTime{1});
    %jumpTimes = time(transition == 1); 
    figure();patch('XData',xPos,'YData',yPos,'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
    hold on
    if s
        patch('XData',xPos(stim == 1),'YData',yPos(stim == 1),'EdgeColor','r','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',7);
    end
    patch('XData',xPos(1),'YData',yPos(1),'EdgeColor','g','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',20);
    %patch(xPos(transition == 1),yPos(transition == 1),time(transition == 1),'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',7);
    set(gcf,'color','w');
    xlabel('mm')
    xlim([min(minX,minY),max(maxX, maxY)])
    ylim([min(minX,minY),max(maxX, maxY)])

    %if savePlot
        %saveas(gcf, fullfile(folder,'figures',['pathTrajectory_trial_',num2str(trial),'.fig']));
        saveas(gcf, fullfile(folder,'figures','pathTrajectoryStim.fig'));
        saveas(gcf, fullfile(folder,'figures','pathTrajectoryStim.svg'));
    %end


end