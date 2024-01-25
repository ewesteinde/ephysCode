%% visualizing heading bump & Vf/Vs relationship
clear

rootPath = '/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys'; % 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'; %change dep on comp
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

if isfield(trialMeta, 'notes')
    disp(trialMeta.notes)
end

cd('/Users/elenawesteinde/Documents/EphysCode/Analysis_code'); 
%%
tStart = 1;
tEnd = 3;
activity_values = [];
x_values = [];
y_values = [];

for type = 1:2
    for t = tStart:tEnd
        
        tData = processed_trialData{t};
    
        if type == 1
            data = getfield(tData,'fRate_sec');
        else
            data = getfield(tData,'smooth_Vm');
        end
        
        activity_values = cat(1,activity_values, data);
        x_values = cat(1, x_values, processed_behaviourData{t}.angle);
        y_values = cat(1, y_values, (processed_behaviourData{t}.vel_for)');

    end

    x_edges = [-180:10:180]; %
    y_edges = [-2:.5:10];


    if type == 1
        [N_fRate, heatmap_fRate, x_centers_fRate, y_centers_fRate] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    else
        [N_Vm, heatmap_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    end
 
end

 figure(3); clf;
            f(1) = subplot(2,2,1);
            imagesc(flip(heatmap_fRate))
            ylabel('Vf mm/s')
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
            
            f(2) = subplot(2,2,3);
            imagesc(flip(heatmap_Vm))
            ylabel('Vf mm/s')
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
            
            hold on
            
            
%
activity_values = [];
x_values = [];
y_values = [];

for type = 1:2
    for t = tStart:tEnd
        
        tData = processed_trialData{t};
    
        if type == 1
            data = getfield(tData,'fRate_sec');
        else
            data = getfield(tData,'smooth_Vm');
        end
        
        activity_values = cat(1,activity_values, data);
        x_values = cat(1, x_values, (processed_behaviourData{t}.vel_side)');
        y_values = cat(1, y_values, (processed_behaviourData{t}.vel_for)');

    end

    x_edges = [-10:.5:10]; 
    y_edges = [-2:1:10];


    if type == 1
        [N_fRate, heatmap_fRate, x_centers_fRate, y_centers_fRate] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    else
        [N_Vm, heatmap_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    end
 
end

        f(3) = subplot(2,2,2);
        imagesc(flip(heatmap_fRate))
        ylabel('Vf mm/s')
        xt = get(gca, 'XTick');
        xtnew = linspace(min(xt), max(xt), 10);                             
        xtlbl = linspace(-8, 10, numel(xtnew));                  
        set(gca, 'XTick',xtnew, 'XTickLabel',xtlbl)
        xlabel('Vs mm/s')
        colorbar
        title('fRate')
        yt = get(gca, 'YTick');
        ytnew = linspace(min(yt),max(yt),6);
        ytlbl = linspace(8.5, -1.5, numel(ytnew));
        set(gca, 'YTick',ytnew, 'YTickLabel',ytlbl) 


        f(4) = subplot(2,2,4);
        imagesc(flip(heatmap_Vm))
        ylabel('Vf mm/s')
        xt = get(gca, 'XTick');
        xtnew = linspace(min(xt), max(xt), 10);                             
        xtlbl = linspace(-8, 10, numel(xtnew));                  
        set(gca, 'XTick',xtnew, 'XTickLabel',xtlbl)
        xlabel('Vs mm/s')
        colorbar
        title('Vm')
        yt = get(gca, 'YTick');
        ytnew = linspace(min(yt),max(yt),6);
        ytlbl = linspace(8.5, -1.5, numel(ytnew));
        set(gca, 'YTick',ytnew, 'YTickLabel',ytlbl) 
%  

% uncomment to check yaw relationship

activity_values = [];
x_values = [];
y_values = [];

for type = 1:2
    for t = tStart:tEnd
        
        tData = processed_trialData{t};
    
        if type == 1
            data = getfield(tData,'fRate_sec');
        else
            data = getfield(tData,'smooth_Vm');
        end
        
        activity_values = cat(1,activity_values, data);
        x_values = cat(1, x_values, (processed_behaviourData{t}.vel_yaw)');
        y_values = cat(1, y_values, (processed_behaviourData{t}.vel_for)');

    end

    x_edges = [-60:5:60]; 
    y_edges = [-2:1:10];


    if type == 1
        [N_fRate, heatmap_fRate, x_centers_fRate, y_centers_fRate] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    else
        [N_Vm, heatmap_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
    end
 
end

figure(2);
        h(1) = subplot(2,1,1);
        imagesc(flip(heatmap_fRate))
        ylabel('Vf mm/s')
%         xt = get(gca, 'XTick');
%         xtnew = linspace(min(xt), max(xt), 10);                             
%         xtlbl = linspace(-8, 10, numel(xtnew));                  
%         set(gca, 'XTick',xtnew, 'XTickLabel',xtlbl)
        xlabel('Vy mm/s')
        colorbar
        title('fRate')
        yt = get(gca, 'YTick');
        ytnew = linspace(min(yt),max(yt),6);
        ytlbl = linspace(8.5, -1.5, numel(ytnew));
        set(gca, 'YTick',ytnew, 'YTickLabel',ytlbl) 


        h(2) = subplot(2,1,2);
        imagesc(flip(heatmap_Vm))
        ylabel('Vf mm/s')
%         xt = get(gca, 'XTick');
%         xtnew = linspace(min(xt), max(xt), 10);                             
%         xtlbl = linspace(-8, 10, numel(xtnew));                  
%         set(gca, 'XTick',xtnew, 'XTickLabel',xtlbl)
        xlabel('Vy mm/s')
        colorbar
        title('Vm')
        yt = get(gca, 'YTick');
        ytnew = linspace(min(yt),max(yt),6);
        ytlbl = linspace(8.5, -1.5, numel(ytnew));
        set(gca, 'YTick',ytnew, 'YTickLabel',ytlbl) 
        
        

keep = input('Save? ','s');

if strcmp(keep, 'y')
    cd(fileName) 
    saveas(gcf,'directionPref.fig')
    cd('/Users/elenawesteinde/Documents/EphysCode/Analysis_code') %'C:\Code\EphysCode'
end


%[-175,-165,-155,-145,-135,-125,-115,-105,-95,-85,-75,-65,-55,-45,-35,-25,-15,-5,5,15,25,35,45,55,65,75,85,95,105,115,125,135,145,155,165,175]
%% If unclear look at trials by eye


if trialMeta.fly.timestamp < '26-Mar-2021' && trialMeta.fly.timestamp > '01-Jan-2021'
    
    cd(fileName)
    load('trialData.mat');
    downsample_Hz = 1000; 
    ephysSettings
    cd('/Users/elenawesteinde/Documents/EphysCode/Analysis_code')
    
    nActivity = resample(trialData{t}.scaledOutput, downsample_Hz, settings.sampRate);
else
    nActivity = processed_trialData{t}.scaledOutput_down;
end

figure(2);clf;      
        
        h(1) = subplot(3,1,1);
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.angle, 'k') 
        ylabel('angle')
       
        
        h(2) = subplot(3,1,2);
        yyaxis left
        plot(processed_behaviourData{t}.time,nActivity, 'k') 
        ylabel('Vm')
        %ylim([-55 -])
        hold on 
        yyaxis right
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_for, 'r')
        ylim([-6 6])
        ylabel('Vf mm/sec')
        
        h(3) = subplot(3,1,3);
        yyaxis left
        plot(processed_behaviourData{t}.time,nActivity, 'k')
        ylabel('Vm')
        %ylim([-55 -])
        hold on 
        yyaxis right
        plot(processed_behaviourData{t}.time, processed_behaviourData{t}.vel_side, 'r')
        ylim([-6 6])
        ylabel('Vs mm/sec')
        
        
        linkaxes(h,'x');