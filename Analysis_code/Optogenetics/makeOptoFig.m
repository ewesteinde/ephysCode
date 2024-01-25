ephysSettings
niIO.Rate = settings.sampRate;
bufferLength = 1;
aveTrial = zeros(length(trialData{1}.scaledOutput),1);
for t = 1:length(trialData)
trialMeta.trialAveBaseline(t) = sum(trialData{t}.scaledOutput(1:bufferLength*niIO.Rate))/length(trialData{t}.scaledOutput(1:bufferLength*niIO.Rate));
trialMeta.trialResp(t) = -(max(abs(trialData{t}.scaledOutput(bufferLength*niIO.Rate:(bufferLength+.2)*niIO.Rate)))-abs(trialMeta.trialAveBaseline(t)));  
aveTrial = aveTrial + (trialData{t}.scaledOutput - trialMeta.trialAveBaseline(t));
end
aveTrial = aveTrial/length(trialData); 
% Plot iclamp trial    
    
aveResp = sum(trialMeta.trialResp(15:end))/length(trialData);
    


figure(11); clf;

                h(1) = subplot(5,1,1:4);
                plot(trialData{t}.time, aveTrial, 'k');
                ylabel('Voltage (mV)')
                ylim([-6 1])
                xlim([0 2])
                
                h(2) = subplot(5,1,5);
                plot(trialData{t}.time, trialData{t}.input,'k');
                ylabel('Light Stim')
                xlabel('Time (s)')
                xlim([0 2])

                sgtitle(['Trial ' num2str(t)])
              
                linkaxes(h,'x')
                
%                 figure(22); clf; 
%                 y = trialMeta.trialAveBaseline(1:t);
%                 yyaxis left
%                 plot(1:t,y,'-.')
%                 hold on
%                 plot(1:t, movmean(trialMeta.trialAveBaseline(1:t),7),'k-');
%                 ylabel('Baseline (mV)')
% 
%                 z = trialMeta.trialResp(1:t);
%                 
%                 yyaxis right
%                 plot(1:t,z,'-.')
%                 ylim([min(trialMeta.trialResp(1:t))-2 0])
%                 hold on
%                 plot(1:t, movmean(trialMeta.trialResp(1:t),7),'r-');
%                 ylabel('Response (mV)')
%                 xlabel('Trial')

%%
figure(5);clf;
sz = 2; 
c = [0 0 1; 0 0 1; 0 0 1; 0 1 0; 0 1 0; 0 1 0];
scatter([1 1 1 1 1 1],aveStimControl(:,3),sz,c,'filled')
ylim([-10 0])
hold on 
scatter([1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5],aveStim(:,4),sz,'filled')
xlim([0.5 2])
ylabel('?Voltage (mV)')
names = {'Control';'Chr Stim'};
set(gca,'xtick',[1,1.5],'xticklabel',names)