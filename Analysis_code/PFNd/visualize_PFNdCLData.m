clear
rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys'; % ;'/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys' 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'%change dep on comp
date = input('Date? ','s');
cell_num = input('Cell? ','s');
cell_num = strcat('cell_',cell_num);
trial = input('Trial? ','s');
trial = strcat('trial_',trial);
fileName = fullfile(rootPath,date,cell_num,trial);

cd(fileName)
load('pro_trialData.mat');
load('pro_behaviourData.mat');
load('trialMeta.mat');
load('trialData.mat');

if isfield(trialMeta, 'notes')
    disp(trialMeta.notes)
end

cd('C:\Code\EphysCode'); %  'C:\Code\EphysCode' '/Users/elenawesteinde/Documents/EphysCode'%change dep on comp


%% Look at entirety of trial & compare VelDot & Vm

    angle = [];
    velDot = [];
    time = [];
    scaledOutput = []; 
count = 0; 
for t = 1:length(processed_trialData) 
    
    angle = cat(1, angle, processed_behaviourData{t}.angle');
    velDot = cat(1, velDot, processed_behaviourData{t}.vel_dot);
    time = cat(1, time, (processed_behaviourData{t}.time+ 300*count)');
    scaledOutput = cat(1,scaledOutput, processed_trialData{t}.scaledOutput_down'); 
    
    count = count + 1; 
end
    
% visualize
% 
figure(2);clf;      
        
        h(1) = subplot(2,1,1);
        plot(time, angle, 'k') 
        ylabel('angle')
       
        h(2) = subplot(2,1,2);
        yyaxis left
        plot(time, scaledOutput, 'k') 
        ylabel('Vm')
        ylim([-70 -35])
        hold on 
        yyaxis right
        plot(time, velDot, 'r')
        %ylim([-6 2])
        
        
        linkaxes(h,'x');
        
        
%% Make rep trace
start = 450;
finish = 467;
figure();clf;
        yyaxis left
       plot(time, scaledOutput, 'k') 
        ylabel('Vm')
        xlim([start finish])
        ylim([-70 -35])
        hold on 
        yyaxis right
        plot(time, velDot, 'r')
        ylabel('velocity (mm/s)')
        xlabel('time (s)')
        ylim([-6 8])
        xlim([start finish])  
        %title('120420 c2 t1_1')
        set(gcf, 'color', 'w')
        
