clear
rootPath = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys'; % ;'/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys' 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'%change dep on comp
 % ;'/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys' 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'%change dep on comp

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

cd('C:\Code\EphysCode'); 

%% separate idx at which fly is moving in CL vs OL

% startOL = 179.368*1000; % sec 070821 cell 3 
% endOL = 204.389*1000;
if iscell(processed_trialData)
    processed_trialData = processed_trialData{1}; 
end

if iscell(processed_behaviourData)
     processed_behaviourData = processed_behaviourData{1}; 
end
%% visualize trial to set startOL time

nActivity = processed_trialData.scaledOutput_down;
figure();clf;      
        
        h(1) = subplot(4,1,1);
        plot(processed_behaviourData.time, processed_behaviourData.angle, 'k') 
        ylabel('angle')
       
        
        h(2) = subplot(4,1,2);
        yyaxis left
        plot(processed_behaviourData.time,nActivity, 'k') 
        ylabel('Vm')
        %ylim([-70 -45])
        hold on 
        yyaxis right
        plot(processed_behaviourData.time, processed_behaviourData.vel_for, 'r')
        ylim([-(max(processed_behaviourData.vel_for)) max(processed_behaviourData.vel_for)])
        ylabel('Vf mm/sec')
        
        
        
        h(3) = subplot(4,1,3);
        yyaxis left
        plot(processed_behaviourData.time,nActivity, 'k')
        ylabel('Vm')
       % ylim([-70 -45])
        hold on 
        yyaxis right
        plot(processed_behaviourData.time, processed_behaviourData.vel_side, 'r')
        ylim([-(max(processed_behaviourData.vel_side)) max(processed_behaviourData.vel_side)])
        ylabel('Vs mm/sec')
        
        h(4) = subplot(4,1,4);
        yyaxis left
        plot(processed_behaviourData.time,nActivity, 'k')
        ylabel('Vm')

        %ylim([-55 -])
        hold on 
        yyaxis right
        plot(processed_behaviourData.time, processed_behaviourData.vel_yaw, 'r')
        ylim([(min(processed_behaviourData.vel_yaw)) max(processed_behaviourData.vel_yaw)])
        ylabel('Vy deg/sec')
        
        
        linkaxes(h,'x');
%%
startOL = 178.188*1000; % sec
endOL = 203.242*1000;
lengthOL =  ((endOL/1000)-(startOL/1000))*1000; 
lengthCL = 179.99*1000;

vf = processed_behaviourData.vel_for;

idx = 1; 
idx_CL = [];
idx_OL = []; 
while idx <= length(vf)
    if idx == 1
        tempCL = 1:startOL; 
        tempOL = max(tempCL)+1:max(tempCL)+lengthOL;
    elseif idx + lengthCL + lengthOL > length(vf)
        if idx + lengthCL >= length(vf)
            tempCL = idx:length(vf);
            tempOL = []; 
        else
            tempCL = idx:idx+lengthCL-1; 
            tempOL = max(tempCL)+1:length(vf);
        end
    else 
        tempCL = idx:idx+lengthCL-1; 
        tempOL = max(tempCL)+1:max(tempCL)+lengthOL; 
    end
    idx = max(max(tempOL), max(tempCL))+1; 
    idx_CL = [idx_CL,tempCL];
    idx_OL = [idx_OL, tempOL];
end

CL.vel_for = processed_behaviourData.vel_for(idx_CL); 
CL.vel_side = processed_behaviourData.vel_side(idx_CL); 
CL.vel_yaw = processed_behaviourData.vel_yaw(idx_CL); 
CL.angle = processed_behaviourData.angle(idx_CL); 

CL.scaledOutput_down = processed_trialData.scaledOutput_down(idx_CL); 
CL.smooth_Vm = processed_trialData.smooth_Vm(idx_CL); 
CL.fRate_sec = processed_trialData.fRate_sec(idx_CL); 

OL.vel_for = processed_behaviourData.vel_for(idx_OL); 
OL.vel_side = processed_behaviourData.vel_side(idx_OL); 
OL.vel_yaw = processed_behaviourData.vel_yaw(idx_OL); 
OL.angle = processed_behaviourData.angle(idx_OL); 

OL.scaledOutput_down = processed_trialData.scaledOutput_down(idx_OL); 
OL.smooth_Vm = processed_trialData.smooth_Vm(idx_OL); 
OL.fRate_sec = processed_trialData.fRate_sec(idx_OL); 

% plot heading vectors of the fly across CL epochs
window = 60; %seconds
minVel = 2; %mm/s
sampRate = 1000; %Hz 

[rho, theta, idx_windows] = PFL2_3_headingDist(window, minVel, CL, sampRate, fileName);

%% split trial data into periods when the fly is likely menotaxing vs not

% currently takes up to 5 min, figure out how to speed up, culprit is when
% I iterate through each meno vs not meno window to parse data into new
% structures


meno_angle = -100; % Ave angle fly seems to be trying to maintain
meno_range = 60; % range around that angle
meno_strength = 0.6; % min heading vector strength to count as menotaxis
not_meno_strength = 0.4;

lim1 = meno_angle - meno_range/2;
lim2 = meno_angle + meno_range/2;

angles = rad2deg(theta); 

if abs(lim1) > 180
    [meno_index] = find(angles >= wrapTo180(lim1) | angles <= lim2 & rho >= meno_strength);
    upperlim = lim2;
    lowerlim = wrapTo180(lim1);
elseif abs(lim2) > 180
    [meno_index] = find(angles >= lim1 | angles <= wrapTo180(lim2) & rho >= meno_strength);
    upperlim = wrapTo180(lim2);
    lowerlim = lim1;
else
    [meno_index] = find(angles >= lim1 & angles <= lim2 & rho >= meno_strength);
    upperlim = lim2;
    lowerlim = lim1;
end

meno_windows = idx_windows(meno_index,:); 

notMeno_index = find(rho < not_meno_strength); 
notMeno_windows = idx_windows(notMeno_index,:);  

menoData.vel_for = [];
menoData.vel_side = [];
menoData.vel_yaw = [];
menoData.angle = [];

menoData.scaledOutput_down = [];
menoData.smooth_Vm = []; 
menoData.fRate_sec = [];

not_menoData.vel_for = [];
not_menoData.vel_side = []; 
not_menoData.vel_yaw = [];
not_menoData.angle = [];

not_menoData.scaledOutput_down = [];
not_menoData.smooth_Vm = [];
not_menoData.fRate_sec = [];

for idx = 1:length(meno_windows(:,1))
    menoData.vel_for = [menoData.vel_for, CL.vel_for(meno_windows(idx,1):meno_windows(idx,2))];
    menoData.vel_side = [menoData.vel_side, CL.vel_side(meno_windows(idx,1):meno_windows(idx,2))];
    menoData.vel_yaw = [menoData.vel_yaw, CL.vel_yaw(meno_windows(idx,1):meno_windows(idx,2))];
    menoData.angle = [menoData.angle, CL.angle(meno_windows(idx,1):meno_windows(idx,2))];

    menoData.scaledOutput_down = [menoData.scaledOutput_down, CL.scaledOutput_down(meno_windows(idx,1):meno_windows(idx,2))];
    menoData.smooth_Vm = [menoData.smooth_Vm, CL.smooth_Vm(meno_windows(idx,1):meno_windows(idx,2))]; 
    menoData.fRate_sec = [menoData.fRate_sec, CL.fRate_sec(meno_windows(idx,1):meno_windows(idx,2))];
end

for idx = 1:length(notMeno_windows(:,1))
    not_menoData.vel_for = [not_menoData.vel_for, CL.vel_for(notMeno_windows(idx,1):notMeno_windows(idx,2))];
    not_menoData.vel_side = [not_menoData.vel_side, CL.vel_side(notMeno_windows(idx,1):notMeno_windows(idx,2))]; 
    not_menoData.vel_yaw = [not_menoData.vel_yaw, CL.vel_yaw(notMeno_windows(idx,1):notMeno_windows(idx,2))];
    not_menoData.angle = [not_menoData.angle, CL.angle(notMeno_windows(idx,1):notMeno_windows(idx,2))];

    not_menoData.scaledOutput_down = [not_menoData.scaledOutput_down, CL.scaledOutput_down(notMeno_windows(idx,1):notMeno_windows(idx,2))];
    not_menoData.smooth_Vm = [not_menoData.smooth_Vm, CL.smooth_Vm(notMeno_windows(idx,1):notMeno_windows(idx,2))];
    not_menoData.fRate_sec = [not_menoData.fRate_sec, CL.fRate_sec(notMeno_windows(idx,1):notMeno_windows(idx,2))];
end

%% compare & contrast analyses of meno vs not meno chunks

iStart_meno = 1; 
iEnd_meno = length(menoData.vel_for);

iStart_notMeno = 1; 
iEnd_notMeno = length(not_menoData.vel_for);


% Activity-behaviour line plots


prefHead = 150; 
range = 120;
step = 0.5;
tStart = 1;
tEnd = 1; 

PFL2_3_behaviourVSactivity_lineplots(tStart, tEnd, prefHead, step, range, date, cell_num, trial, menoData, menoData, fileName, iStart_meno, iEnd_meno) ;
PFL2_3_behaviourVSactivity_lineplots(tStart, tEnd, prefHead, step, range, date, cell_num, trial, not_menoData, not_menoData, fileName, iStart_notMeno, iEnd_notMeno) ;

figure();
subplot
histogram(menoData.smooth_Vm,'BinWidth',0.2);
hold on; 
histogram(not_menoData.smooth_Vm,'BinWidth',0.2)

%%

iStart_meno = 1; 
iEnd_meno = 6000;

iStart_notMeno = 30001; 
iEnd_notMeno = 36000;

basic_xcorr_PFL(trialData, menoData, menoData, date, cell_num, trial, fileName,iStart_meno,iEnd_meno) 
basic_xcorr_PFL(trialData, not_menoData, not_menoData, date, cell_num, trial, fileName,iStart_notMeno,iEnd_notMeno)

%% look at heading tuning in OL chunks following high meno & low meno CL epochs

% doesn't work b/c doesn't account for the fact that OL segments happen
% every 3 min not every one, just do it manually for now
% OL_Meno_idx = [];
% OL_noMeno_idx = [];
% 
% for min = 1:length(notMeno_index)
%     startOL_temp = startOL + lengthOL*(notMeno_index(min)-1) + lengthCL*(notMeno_index(min)-1);
%     endOL_temp = startOL + lengthOL*notMeno_index(min) + lengthCL*(notMeno_index(min)-1);
%     OL_noMeno_idx = [OL_noMeno_idx, [startOL_temp:endOL_temp]];
% end
% 
% for min = 1:length(meno_index)
%     startOL_temp = startOL + lengthOL*(meno_index(min)-1) + lengthCL*(meno_index(min)-1);
%     endOL_temp = startOL + lengthOL*meno_index(min) + lengthCL*(meno_index(min)-1);
%     OL_Meno_idx = [OL_Meno_idx, [startOL_temp:endOL_temp]];
% end

% 070121 C1 T1 OL indices 
OL_idx{1} = [startOL:endOL]; % follows CL min 1-3
OL_idx{2} = [383.222*1000:408.255*1000]; % follows CL min 4-6
OL_idx{3} = [588.25*1000:613.246*1000]; % follows CL min 7-9
OL_idx{4} = [793.255*1000:818.26*1000]; % follows CL min 9-12
OL_idx{5} = [998.236*1000:1023.27*1000]; % follows CL min 13-15

%  startOL = 179.368*1000; % sec 070821 cell 3 
%  endOL = 204.389*1000;
% 
% OL_idx{1} = [startOL:endOL]; % follows CL min 1-3
% OL_idx{2} = [384.401*1000:409.395*1000]; % follows CL min 4-6
% OL_idx{3} = [589.382*1000:614.424*1000]; % follows CL min 7-9
% OL_idx{4} = [794.426*1000:819.403*1000]; % follows CL min 10-12
% OL_idx{5} = [999.405*1000:1024.421*1000]; % follows CL min 13-15
% OL_idx{6} = [1204.42*1000:1229.43*1000]; % follows CL min 16-18
% OL_idx{7} = [1409.5*1000:1434.44*1000]; % follows CL min 19-21
% OL_idx{8} = [1614.4*1000:1639.48*1000]; % follows CL min 22-24


for OL_epoch = 1:length(OL_idx)
    OL_temp.smooth_Vm = processed_trialData.smooth_Vm(OL_idx{OL_epoch});
    OL_temp.fRate_sec = processed_trialData.fRate_sec(OL_idx{OL_epoch});
    OL_temp.angle = processed_behaviourData.angle(OL_idx{OL_epoch});
    OL_temp.vel_for = processed_behaviourData.vel_for(OL_idx{OL_epoch});
    OL_temp.speed = sqrt(abs(processed_behaviourData.vel_for(OL_idx{OL_epoch})).^2 + abs(processed_behaviourData.vel_side(OL_idx{OL_epoch})).^2);
    
    figure();clf;
    Vm = OL_temp.smooth_Vm;
    angle = OL_temp.angle;
    edges = [-180:30:180];
    [centers, mean_bin] = create_binned_mean(angle, Vm, edges);
    plot(centers, mean_bin,'-o');
    xlabel('Angle');
    ylabel('Vm'); 
    title('OL heading pref')
    
    
    for a = 1:2
        activity_values = [];
        x_values = [];
        y_values = [];
        for t = tStart:tEnd

            tData = OL_temp;

            if a == 1
                data = getfield(tData,'fRate_sec');
            else
                data = getfield(tData,'smooth_Vm');
            end

            activity_values = cat(1,activity_values, data');
            x_values = cat(1, x_values, OL_temp.angle');
            y_values = cat(1, y_values, (OL_temp.speed'));

        end

        x_edges = [-180:30:180]; %
        y_edges = [-2:.5:10];


        if a == 1
            [N_fRate, heatmap_fRate, x_centers_fRate, y_centers_fRate] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
        else
            [N_Vm, heatmap_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
        end

    end

 figure(); clf;
            f(1) = subplot(2,1,1);
            imagesc(flip(heatmap_fRate))
            ylabel('Speed mm/s')
            xt = get(gca, 'XTick');
            xtnew = linspace(min(xt), max(xt), 7);                             
            xtlbl = linspace(-135, 165, numel(xtnew));                  
            set(gca, 'XTick',xtnew, 'XTickLabel',xtlbl)
            xlabel('angle')
            colorbar
            title('fRate')
            yt = get(gca, 'YTick');
            ytnew = linspace(min(yt),max(yt),6);
            ytlbl = linspace(8.5, -1.5, numel(ytnew));
            set(gca, 'YTick',ytnew, 'YTickLabel',ytlbl)
            
            f(2) = subplot(2,1,2);
            imagesc(flip(heatmap_Vm))
            ylabel('Speed mm/s')
            xt = get(gca, 'XTick');
            xtnew = linspace(min(xt), max(xt), 7);                             
            xtlbl = linspace(-135, 165, numel(xtnew));                  
            set(gca, 'XTick',xtnew, 'XTickLabel',xtlbl)
            xlabel('angle')
            colorbar
            title('Vm')
            yt = get(gca, 'YTick');
            ytnew = linspace(min(yt),max(yt),6);
            ytlbl = linspace(8.5, -1.5, numel(ytnew));
            set(gca, 'YTick',ytnew, 'YTickLabel',ytlbl)
            set(gcf,'color','w')
            
end
