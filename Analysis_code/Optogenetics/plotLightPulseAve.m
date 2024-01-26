%% Plot ave stim response over all trials
trialDataAve = zeros(size(trialData{1, 1}.scaledOutput));
for t = 1:length(trialData)
    trialDataAve = trialDataAve + trialData{1,t}.scaledOutput;
end
trialDataAve = trialDataAve/length(trialData);



figure; clf;
h(1) = subplot(4,1,1:3);
plot(trialData{t}.time, trialDataAve, 'k')
ylabel('Voltage (mV)')

h(2) = subplot(4,1,4);
plot(trialData{t}.time, trialData{t}.input,'k')
ylabel('Light Stim')
xlabel('Time (s)')
sgtitle(['TrialAve ' num2str(trial)])
linkaxes(h,'x')

%% Plot stim resp/baseline trend 
figureNum = 3;
t=trialMeta.trials;
timeBase = [0:4:4*length(totalBase)-1]./60;
timeStim = [2:4:4*length(totalBase)+1]./60;

figure(figureNum); clf;
y = totalBase;%trialMeta.trialAveBaseline(1:t);
yyaxis left
plot(timeBase,y,'-.')
hold on
plot(timeBase, movmean(totalBase,7),'k-');
ylim([-40 -20])
ylabel('Baseline (mV)')

z = totalstim;%trialMeta.trialResp(1:t);

yyaxis right
plot(timeStim,z,'-.')
ylim([min(trialMeta.trialResp(1:t))-2 0])
hold on
plot(timeStim, movmean(totalstim,7),'r-');
xline(0,'-',{'TTX'})
xline(596/60,'-',{'1uM Picrotoxin'})
xline(1156/60,'-',{'100uM Picrotoxin'})
ylabel('Response (mV)')
xlabel('Time(min)')
%title('101920 1uM TTX picro contd')

%% Plot stim ave trend across exps
nCells = 5;
nStates = 2;
figData = cell(nStates, nCells);
for c = 1:length(figData(1,:))
    for s = 1:length(figData(:,1))
        trial = input('What cell & trial folder? ','s');
        cd(trial) 
        fileName = strcat(trial,'\trialMeta.mat');
        load(fileName)
        cd C:\Code\JennyCode
        figData(s,c) = {trialMeta.trialResp};
    end
end


aveResp = zeros(size(figData{1}));
for c = 1:length(figData(1,:))
    aveResp = aveResp + figData{2,c};
end
aveResp = aveResp/nCells;

figure(2); clf; 
y(1) = 2;
y(2:150) = 2+4*(1:(length(figData{1})-1));
%plot(y,figData{2,1},'-.') 
plot(y, movmean(figData{2,1},7),'.');
hold on
%plot(y,aveResp)
plot(y,movmean(aveResp,7),'k','LineWidth',2)
for c = 2:length(figData(1,:))
    %plot(y,figData{2,c},'-.')
    plot(y, movmean(figData{2,c},7),'.');
end

ylim([-10 0])
ylabel('Response (mV)')
xlabel('Time (s)')
title('test')
%%
aveResp = zeros(size(cpgfigData{1}));
for c = 1:length(cpgfigData(1,:))
    aveResp = aveResp + cpgfigData{1,c};
end
aveResp = aveResp/3;

figure(4); clf; 
y(1) = 2;
y(2:150) = 2+4*(1:(length(cpgfigData{1})-1));
%plot(y,figData{2,1},'-.') 
plot(y, movmean(cpgfigData{1,1},7),'.');
hold on
%plot(y,aveResp)
plot(y,movmean(aveResp,7),'k','LineWidth',2)
for c = 2:length(cpgfigData(1,:))
    %plot(y,figData{2,c},'-.')
    plot(y, movmean(cpgfigData{1,c},7),'.');
end

ylim([-10 0])
ylabel('Response (mV)')
xlabel('Time (s)')
title('test')
%%
catfigData = cell(1,nCells);
for c = 1:nCells
    catfigData(:,c) = {cat(2,figData{1,c},figData{2,c},figData{3,c})};
end 
empty = ones(1,150);
aveResp = zeros(size(catfigData{1}));
for c = 1:length(catfigData(1,:))
    x = catfigData{1,c};
    aveResp(1:450) = aveResp(1:450) + x(1:450);
end
aveResp(1:450) = aveResp(1:450)/5;
for c = 2:length(catfigData(1,:))-1
    b = catfigData{1,c};
    aveResp(451:600) = aveResp(451:600) + b(451:600);
end
aveResp(451:600) = aveResp(451:600)/3;

figure(5); clf;
y2 = 1:length(catfigData{1});
plot(y2,movmean(catfigData{1},7),'.')
hold on 
plot(y2,movmean(aveResp,7),'k','LineWidth',2)
for c = 2:length(catfigData(1,:))
    %plot(y,figData{2,c},'-.')
    plot(y2, movmean(catfigData{1,c},7),'.');
end

ylim([-10 0])
ylabel('Response (mV)')
xlabel('Trial')
title('test')
