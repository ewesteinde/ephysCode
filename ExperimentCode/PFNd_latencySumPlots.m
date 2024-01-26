FRlags = [0.0047, -0.2073, -0.0443, -0.0803, -0.0413, 0.0363, -0.0470, -3.3333e-04, 0.0263, 0.0350, 0.0243, 0.0350, 0.0323];
ave_FRlags = mean(FRlags);
med_FRlags = median(FRlags);

Vmlags = [0.0077, -0.0087, -0.1253, -0.0700, -0.0940, 0.0530,  -0.0353, 0.0043, 0.0383, 0.1670, 0.0567, 0.0277, 0.0707, 0.0447];
ave_Vmlags = mean(Vmlags);
med_Vmlags = median(Vmlags);

FR_no0Vel_lags = [-0.0483, -0.0813, -0.0630, -0.0073, -0.0690, -0.0300, 0.0033, 0.0127, -0.0183, 0.0090, -0.0013];
ave_FR_no0Vel_lags = mean(FR_no0Vel_lags);
med_FR_no0Vel_lags = median(FR_no0Vel_lags);

y_FRlags = ones(length(FRlags)); 
y_Vmlags = ones(length(Vmlags)); 
y_FR_no0Vel_lags = ones(length(FR_no0Vel_lags)); 

figure(1); clf; 
plot(FRlags, y_FRlags,'ko')
hold on
plot(med_FRlags, 1, 'ro','MarkerFaceColor','r')
xlim([-0.25 0.25])
ylim([0.5 1.5])
xlabel('Time of peak speed-firing rate correlation (s)')
set(gca, 'YTick', [])
set(gcf,'color','w')

figure(2); clf; 
plot(Vmlags, y_Vmlags,'ko')
hold on
plot(med_FRlags, 1, 'ro','MarkerFaceColor','r')
xlim([-0.25 0.25])
ylim([0.5 1.5])
xlabel('time of peak speed-voltage correlation (s)')
set(gca, 'YTick', [])
set(gcf,'color','w')

figure(3); clf; 
plot(FR_no0Vel_lags, y_FR_no0Vel_lags,'ko')
hold on 
plot(med_FR_no0Vel_lags, 1, 'ro','MarkerFaceColor','r')
xlim([-0.25 0.25])
ylim([0.5 1.5])
xlabel('Time of peak speed-firing rate correlation within epochs of continuous movement (s)')
set(gca, 'YTick', [])
set(gcf,'color','w')

