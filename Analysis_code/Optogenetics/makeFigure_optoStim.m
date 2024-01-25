%% plot ave of trials
start = 1;
last1 = length(trialData);
trialDataAve1 = zeros(size(trialData{1, 1}.scaledOutput));
figure; clf;
h(1) = subplot(4,1,1:3);
for t = start:last1
    trialDataAve1 = trialDataAve1 + trialData{t}.scaledOutput;
    plot(trialData{t}.time, trialData{t}.scaledOutput)
    hold on
end
trialDataAve1 = trialDataAve1/t;
plot(trialData{t}.time, trialDataAve1, 'k', 'LineWidth',1.5) 
ylabel('Voltage (mV)')

h(2) = subplot(4,1,4);
plot(trialData{t}.time, trialData{t}.input,'k')
ylabel('Light Stim')
xlabel('Time (s)')
sgtitle(strcat('092420 c4-t2-t',num2str(start),'-',num2str(last1), '-1uM TTX'))
linkaxes(h,'x')

%% plot one trial 
t= 50;

figure; clf;
h(1) = subplot(4,1,1:3);
plot(trialData{t}.time, trialData{t}.voltage, 'k') 
ylabel('Voltage (mV)')

h(2) = subplot(4,1,4);
plot(trialData{t}.time, trialData{t}.input,'k')
ylabel('Light Stim')
xlabel('Time (s)')
sgtitle(strcat('092420 c4-t2-t',num2str(t), '-1uM TTX'))
linkaxes(h,'x')
%% Plot two overlayed ave
start1 = 1;
last1 =60 ;%length(trialData);


trialDataAve1 = zeros(size(trialData{1, 1}.scaledOutput));
for t = start1:last1
    trialDataAve1 = trialDataAve1 + trialData{1,t}.scaledOutput;
end
trialDataAve1 = trialDataAve1/(last1-start1);
trialDataAve1 = trialDataAve1 - mean(trialDataAve1(1:2*(2e4),1));
%%
start2 = 100;
last2 = length(trialData);
trialDataAve2 = zeros(size(trialData{1, 1}.scaledOutput));
for t = start2:last2
    trialDataAve2 = trialDataAve2 + trialData{1,t}.scaledOutput;
end
trialDataAve2 = trialDataAve2/(last2-start2);
trialDataAve2 = trialDataAve2 - mean(trialDataAve2(1:2*(2e4),1));


%%
figure; clf;
%h(1) = subplot(4,1,1:3);
plot(trialData{t}.time(20000:80000), trialDataAve1(20000:80000), 'Color', [0.8 0 0])
hold on
plot(trialData{t}.time(20000:80000), trialDataAve2(20000:80000),'k' )
xlim([1.5 4])
ylim([-8 1])
 x = [2 2.01 2.01 2];
 y = [.8 .8 .5 .5];
% col[1,1,3] = [0.8500 0.3250 0.0980];
patch(x,y,[0.8500 0.3250 0.0980],'EdgeColor','none','faceAlpha',1)
ylabel('?Voltage (mV)')
% legend('ATR','No ATR')
% legend boxoff
box off
% h(2) = subplot(4,1,4);
% plot(trialData{t}.time(20000:80000), trialData{t}.input(20000:80000),'k')
% ylabel('Cy5 Stim')
xlabel('Time (s)')
% %sgtitle(strcat('PFNd Recording SPSP CsChrimson Stim'))
% xlim([1.5 4])
% linkaxes(h,'x')

%% Plot three overlayed ave
start1 = 1;
last1 = 65;
start2 = 85;
last2 = length(trialData);
%%
start3 = 60;
last3 = length(trialData);

trialDataAve3 = zeros(size(trialData{1, 1}.scaledOutput));
for t = start3:last3
    trialDataAve3 = trialDataAve3 + trialData{1,t}.scaledOutput;
end
trialDataAve3 = trialDataAve3/(last3-start3);
trialDataAve3 = trialDataAve3 - mean(trialDataAve3(1:2*(2e4),1));
%%
figure; clf;
h(1) = subplot(4,1,1:3);
plot(trialData{t}.time, trialDataAve1, 'k')
hold on
plot(trialData{t}.time, trialDataAve2, 'r')
plot(trialData{t}.time, trialDataAve3, 'b')
ylabel('?Voltage (mV)')
legend('1 uM TTx', '+5uM picrotoxin','+25uM picro')

h(2) = subplot(4,1,4);
plot(trialData{t}.time, trialData{t}.input,'k')
ylabel('Light Stim')
xlabel('Time (s)')
sgtitle(strcat('092420 Cell4'))
linkaxes(h,'x')