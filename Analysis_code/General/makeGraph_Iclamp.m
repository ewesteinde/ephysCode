fileName = input('what file?','s');
fileName = fullfile(fileName,'trialData.mat');
load(fileName);
k = strfind(fileName,'cell');
cell = fileName(k:k+5);
cell(cell == '_') = ' ';
t = 1;

trialTime =(1:length(trialData{1}(:,1)))/10^5;

% @ sample rate of 10^5 limit of 7 min chunks (420s) saved in a single
% figure as .fig
lengthTrial = max(trialTime)
part = input('Plot the whole trial? y/n ','s');
if ~strcmp(part,'y')
    chunk = input('chunk to plot (s) ');
    trialChunk = trialTime(trialTime>min(chunk) & trialTime<=max(chunk));
    start = min(chunk)*10^5;
    finish = max(chunk)*10^5-1;
else
    trialChunk = trialTime;
    start = 1;
    finish = length(trialData{1}(:,1)); 
end

 % Plot iclamp trial
figure(1); clf;
h(1) = subplot(4,1,1);
plot(trialChunk', trialData{t}(start:finish,1), 'k')
ylabel('Current (pA)')

h(2) = subplot(4,1,2:4);
plot(trialChunk', trialData{t}(start:finish,3), 'k')
ylabel('Voltage (mV)')
xlabel('Time (s)')

sgtitle([cell ' trial ' num2str(t)])

linkaxes(h,'x')

keep = input('Save figure? y/n ','s');

if strcmp(keep, 'y')
    figName = fileName(1:end-13);
    label = input('Figure name? ','s');
    figName = strcat(figName,label,'.fig');
    savefig(figName);
end
%% Load meta files
fileName = input('what file?','s');
fileName = fullfile(fileName,'trialMeta.mat');
load(fileName);
k = strfind(fileName,'cell');
cellNum = fileName(k+5);
resArray(t,1) = str2double(cellNum);
resArray(2,2) = trialMeta.sealR;
resArray(t,3) = trialMeta.initialInputR;
