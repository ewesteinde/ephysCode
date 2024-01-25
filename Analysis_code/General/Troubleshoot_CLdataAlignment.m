%test different smoothing functions

%% check raw data alignment 
  ephysSettings  
% angle is calculated from rawData( :, settings.raw.panelsX )
    %this is sent directly from the controller to the DAQ where its fed
    %into the computer

    raw_angle = rawData{1}( :, settings.raw.panelsX );
% yaw displacement is calculated from rawData( :, settings.raw.yaw)
    % this is sent from the computer to the phidget where its converted
    % back into analog signal & fed to the DAQ
    raw_yaw = rawData{1}( :, settings.raw.yaw);
% if there is no distinct lag between these two traces then the cue
% position is correctly tracking the ball yaw position & any delay in the
% processed traces is coming from the processing itself but the raw data is
% fine
% 03/17/20 cell trace: I see at most a ~20ms delay b/w the raw angle & yaw trace, they follow
% eachother quite well 
figure(1); clf; 
        yyaxis left
       % h(1) = subplot(3,1,1);
        plot(raw_angle);
        ylim([0 10])
        yyaxis right
        %h(2) = subplot(3,1,2);
        plot( raw_yaw);
        ylim([0 10])
 
%%
[ vel_yaw_noSmooth , disp_yaw_noSmooth ] = ficTracSignalDecoding_noSmooth( rawData{1}( :, settings.raw.yaw) , settings.sampRate , settings.fictracRate, 1);
%% smooth displacement before calculating vel
[ vel_yaw_dispSmooth, disp_yaw_dispSmooth ] = ficTracSignalDecoding(rawData{1}( :, settings.raw.yaw) , settings.sampRate , settings.fictracRate, 1);
%based on the signal analyzer smoothdata rloess seemed to adequatly smooth
%the signal with the least temporal spreading 

figure(2); clf; 
        h(1) = subplot(2,1,1);
        plot( vel_yaw_noSmooth,'k');
        hold on
        plot(vel_yaw_dispSmooth,'-r');
        h(2) = subplot(2,1,2);
        plot( disp_yaw_noSmooth,'k');
        hold on
        plot(disp_yaw_dispSmooth,'-r');
        
        linkaxes(h,'x')
%% smooth vel but not displacement 
[ vel_yaw_velSmooth, disp_yaw_velSmooth ] = ficTracSignalDecoding_velSmooth(rawData{1}( :, settings.raw.yaw) , settings.sampRate , settings.fictracRate, 1);

figure(3); clf; 
        h(3) = subplot(2,1,1);
        plot( vel_yaw_noSmooth,'k');
        hold on
        plot(vel_yaw_velSmooth,'-r');
        h(4) = subplot(2,1,2);
        plot( disp_yaw_noSmooth,'k');
        hold on
        plot(disp_yaw_velSmooth,'-r');
        
        linkaxes(h,'x')
        
%% compare gradient to rdiff & fidelity to angle & displacement 

vel_yaw_rdiff = rdiff(disp_yaw_noSmooth, 1/30); %better signal to noise than gradient, takes 6 min per 5 min trace though 
[ ~, ~, angle] = process_panel_360( [rawData{1}( :, settings.raw.panelsX) rawData{1}( :, settings.raw.panelsY )],settings.panels.numFramesX,settings.panels.numFramesY,8, 1);

 fHz = 1000;
    %[y, t] = resample_new(x, fs_new, fs_old)
    [int_vel_yaw_noSmooth,time] = resample_new(vel_yaw_noSmooth,fHz,(settings.fictracRate/2));
    [int_vel_yaw_rdiff,~] = resample_new(vel_yaw_rdiff,fHz,(settings.fictracRate/2));
   
    [int_angle,~] = resample_new(angle,fHz,(settings.sampRate));
    int_angle = smoothdata(int_angle,'movmedian',20);
    int_angle(angle > 180) = 180;
    int_angle(angle < -180) = -180; 


figure(4); clf; 
        j(1) = subplot(2,1,1);
        plot(time, int_vel_yaw_noSmooth,'k');
        hold on
        plot(time,int_vel_yaw_rdiff,'-r');
        j(2) = subplot(2,1,2);
        plot(time,int_angle)
        
        linkaxes(j,'x')

vel_yaw_rdiff_smooth =  smoothdata(vel_yaw_rdiff,'loess',10); %worthwile

[int_vel_yaw_velSmooth,time] = resample_new(vel_yaw_velSmooth,fHz,(settings.fictracRate/2));
[int_vel_yaw_rdiff_smooth,~] = resample_new(vel_yaw_rdiff_smooth,fHz,(settings.fictracRate/2));        
        
figure(5); clf; 
    j(1) = subplot(2,1,1);
    plot(time, int_vel_yaw_velSmooth,'k');
    hold on
    plot(time,int_vel_yaw_rdiff_smooth,'-r');
    j(2) = subplot(2,1,2);
    plot(time,int_angle)

    linkaxes(j,'x')
