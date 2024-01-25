% Example velocity traces

clear

[processed_trialData, processed_behaviourData, fileName] = loadPFL2_3_CLOL('PFL3');

%%

figure(1);clf;      
       
        h(1) = subplot(3,1,1);
        plot(processed_behaviourData.time, processed_behaviourData.angle, 'k') 
        ylabel('angle')
        ylim([-180 180])
        box off
               
        h(2) = subplot(3,1,2);
        yyaxis right
        plot(processed_behaviourData.time, processed_behaviourData.vel_for) 
        yyaxis left 
        plot(processed_behaviourData.time, processed_behaviourData.vel_yaw) 
        ylabel('Speed (mm/s)')
        %ylim([-8 8])
        box off
       
        
        h(3) = subplot(3,1,3);
        plot(processed_behaviourData.time,processed_trialData.scaledOutput_down, 'k') 
        ylabel('Vm (mV)')
        set(gcf,'color','w')
        box off
 
        linkaxes(h,'x');


%% velocity plots

start = 30;
finish = 50;%processed_behaviourData{1}.time(end);
t = 1;

% basic activity vs speed plots
speed_wholeTrial = sqrt(processed_behaviourData.vel_for.^2 + processed_behaviourData.vel_side.^2);

iStart = find(processed_behaviourData.time == start); 
iEnd = find(processed_behaviourData.time == finish); 

nActivity = processed_trialData.scaledOutput_down(iStart:iEnd); 
speed = speed_wholeTrial(iStart:iEnd); 
vf = processed_behaviourData.vel_for(iStart:iEnd); 
vy = processed_behaviourData.vel_yaw(iStart:iEnd); 
time = processed_behaviourData.time(iStart:iEnd); 
angle = processed_behaviourData.angle(iStart:iEnd); 

figure();clf;      
        
        h(1) = subplot(3,1,1);
        plot(time, vf, 'k') 
        ylabel('Forward Velocity (mm/s)')
        ylim([-1 12])
        xlim([start finish])
        set(gca,'XColor','none')
        box off

        h(2) = subplot(3,1,2);
        plot(time,nActivity, 'k') 
        ylabel('Vm (mV)')
        set(gcf,'color','w')
        ylim([-70 -50])
        xlim([start finish])
        set(gca,'XColor','none')
        box off     
        
        h(3) = subplot(3,1,3);
        plot(time, vy, 'k') 
        ylabel('Yaw Velocity (deg/s)')
        ylim([-200 200])
        xlim([start finish])
        y = -100; 
        l = line([finish-2 finish], [y y]);
        l.Color = 'k';
        set(gca,'XColor','none')
        box off

 
        linkaxes(h,'x');
        
%% angle plots

%% velocity plots

start = 626;
finish = 649;%processed_behaviourData{1}.time(end);

% basic activity vs speed plots

iStart = find(processed_behaviourData.time == start); 
iEnd = find(processed_behaviourData.time == finish); 

nActivity = processed_trialData.scaledOutput_down(iStart:iEnd); 
time = processed_behaviourData.time(iStart:iEnd); 
angle = processed_behaviourData.angle(iStart:iEnd); 
vf = processed_behaviourData.vel_for(iStart:iEnd); 
vy = processed_behaviourData.vel_yaw(iStart:iEnd); 

%% alternative overlay plots

fig = figure();clf;
set(gcf,'renderer','painters')
    yyaxis left
    plot(time, nActivity, 'k')
    ylabel('Vm (mV)')
    ylim([-70 -40])
    y = -70;
    l = line([finish-3 finish-1], [y y]);
    l.Color = 'k';
    yyaxis right
    plot(time, angle,'r')
    ylabel('Forward Velocity (mm/sec)')
%      ylim([-180 180])
    xlim([start finish])
    left_color = [0 0 0];
    right_color = [1 0 0];
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
    set(gca,'XColor','none')
    set(gcf,'color',[1 1 1])
    ylim([-200 200])

    box off
%%
figure();clf;      
        
        h(1) = subplot(4,1,1);
        plot(time, angle, 'k') 
        ylabel('Cue Position (deg)')
        ylim([-100 100])
        xlim([start finish])
        set(gca,'XColor','none')
        box off

        h(2) = subplot(4,1,2);
        plot(time,nActivity, 'k') 
        ylabel('Vm (mV)')
        set(gcf,'color','w')
        set(gca,'XColor','none')
        ylim([-70 -45])
        xlim([start finish])
        box off
        
        
         h(3) = subplot(4,1,3);
        plot(time,vf, 'k') 
        ylabel('Forward Vel (mm/s)')
        set(gca,'XColor','none')
        set(gcf,'color','w')
        box off
        
         h(4) = subplot(4,1,4);
        plot(time,vy, 'k') 
        ylabel('Yaw Vel (deg/s)')
        set(gcf,'color','w')
        y = -100; 
        l = line([finish-2 finish], [y y]);
        l.Color = 'k';
        set(gca,'XColor','none')
        box off     
        
 
        linkaxes(h,'x');
        

