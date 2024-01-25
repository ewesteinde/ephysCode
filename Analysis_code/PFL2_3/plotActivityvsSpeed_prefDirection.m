%close all
% clear
% rootPath = 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys' ; % ;  '/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys'%change dep on comp
% date = input('Date? ','s');
% cell_num = input('Cell? ','s');
% cell_num = strcat('cell_',cell_num);
% trial = input('Trial? ','s');
% trial = strcat('trial_',trial);
% fileName = fullfile(rootPath,date,cell_num,trial);
% 
% cd(fileName)
% load('pro_trialData.mat');
% load('pro_behaviourData.mat');
% openfig('directionPref_velDot.fig');
% 
% 
% cd('C:\Code\EphysCode'); %'/Users/elenawesteinde/Documents/EphysCode'

%%
function [saveBinsFRVel, saveBinsVmVel] = plotActivityvsSpeed_prefDirection(tStart, tEnd, prefHead, step, range, offset, maxValPlot, minValPlot, date, cell_num, trial, processed_trialData, processed_behaviourData, fileName) 
%Speed 
    name = strcat(date,', ',cell_num,', ',trial,' prefAngle: ',num2str(prefHead),' BinSize: ',num2str(range),' Offset: ', num2str(offset));
    angles = {num2str(wrapTo180(prefHead+180)), strcat(num2str(wrapTo180(prefHead+90)),'/',num2str(wrapTo180(prefHead-90))), num2str(prefHead)};
    headings = [wrapTo180(prefHead+180),wrapTo180(prefHead+90), prefHead];
    
    
    Vm = [];
    fRate_sec = [];
    angle = [];
    velDot = [];

    for t = tStart:tEnd 

        fRate_sec = cat(1,fRate_sec, processed_trialData{t}.fRate_sec);
        Vm = cat(1,Vm, processed_trialData{t}.smooth_Vm);
        angle = cat(1, angle, processed_behaviourData{t}.angle);
        velDot = cat(1, velDot, processed_behaviourData{t}.vel_dot);
    end
    %figure; histogram(angle)
    
    speedDot = abs(velDot);
    
    count = 1; 
    saveBinsVmS = cell(1,3);
    saveBinsFRS = cell(1,3); 
    
    saveBinsVmVel = cell(1,3);
    saveBinsFRVel = cell(1,3); 
    
    for head = headings
        VinHead = []; 
        VminHead = []; 
        fRateinHead = []; 
        index = [];

        if count == 2
            rangeh = range/2;
            
            lim1 = head-rangeh/2;
            lim2 = head+rangeh/2;
            
            lim3 = wrapTo180(head+180)-rangeh/2;
            lim4 = wrapTo180(head+180)+rangeh/2;
        
            if abs(lim1) > 180
                [index1] = find(angle >= wrapTo180(lim1) | angle <= lim2);
                upperlim1 = lim2;
                lowerlim1 = wrapTo180(lim1);
            elseif abs(lim2) > 180
                [index1] = find(angle >= lim1 | angle <= wrapTo180(lim2));
                upperlim1 = wrapTo180(lim2);
                lowerlim1 = lim1;
            else
                [index1] = find(angle >= lim1 & angle <= lim2);
                upperlim1 = lim2;
                lowerlim1 = lim1;
            end
            
            if abs(lim3) > 180
                [index2] = find(angle >= wrapTo180(lim3) | angle <= lim4);
                upperlim2 = lim4;
                lowerlim2 = wrapTo180(lim3);
            elseif abs(lim4) > 180
                [index2] = find(angle >= lim3 | angle <= wrapTo180(lim4));
                upperlim2 = wrapTo180(lim4);
                lowerlim2 = lim3;
            else
                [index2] = find(angle >= lim3 & angle <= lim4);
                upperlim2 = lim4;
                lowerlim2 = lim3;
            end
            
            index = cat(1,index1, index2); 
        else

            lim1 = head-range/2;
            lim2 = head+range/2;
        
            if abs(lim1) > 180
                [index] = find(angle >= wrapTo180(lim1) | angle <= lim2);
                upperlim = lim2;
                lowerlim = wrapTo180(lim1);
            elseif abs(lim2) > 180
                [index] = find(angle >= lim1 | angle <= wrapTo180(lim2));
                upperlim = wrapTo180(lim2);
                lowerlim = lim1;
            else
                [index] = find(angle >= lim1 & angle <= lim2);
                upperlim = lim2;
                lowerlim = lim1;
            end
        end
        
        VinHead = velDot(index); 
        SinHead = speedDot(index);
        VminHead = Vm(index); 
        fRateinHead = fRate_sec(index);

                edges = [-(step/2):step:max(SinHead(:,1))]; %start at -step/2 so center of first bin is 0mm/s
                [N, edges, bin] = histcounts(SinHead(:,1), edges);
                tempFR = accumarray(bin+1, fRateinHead, [length(edges) 1]);
                tempVm = accumarray(bin+1, VminHead, [length(edges) 1]);
                mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
                mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
                centers = edges(1:end-1)+diff(edges)/2;
                
                saveBinsVmS{count} = [centers' mean_binVm];
                saveBinsFRS{count} = [centers' mean_binFR];
                
                
                edges = [min(VinHead(:,1)):step:max(VinHead(:,1))]; %
                [N, edges, bin] = histcounts(VinHead(:,1), edges);
                tempFR = accumarray(bin+1, fRateinHead, [length(edges) 1]);
                tempVm = accumarray(bin+1, VminHead, [length(edges) 1]);
                mean_binFR = bsxfun(@rdivide, tempFR(2:end), N');
                mean_binVm = bsxfun(@rdivide, tempVm(2:end), N');
                centers = edges(1:end-1)+diff(edges)/2;
                
                saveBinsVmVel{count} = [centers' mean_binVm];
                saveBinsFRVel{count} = [centers' mean_binFR];
        count = count + 1; 
    end
    
% uncomment to plot speed (will need to edit to put in lines of best fit)
% 
%    g(1) = figure(2); clf; 
%     subplot(2,1,1);
%     plot(saveBinsVmS{1}(:,1),saveBinsVmS{1}(:,2))
%     hold on 
%     xlim([0 maxValPlot])
%     plot(saveBinsVmS{2}(:,1),saveBinsVmS{2}(:,2))   
%     plot(saveBinsVmS{3}(:,1),saveBinsVmS{3}(:,2))
%     legend(angles)
%     xlabel('Speed mm/sec')
%     ylabel('Vm (mV)')
%     title(strcat(name, ' Speed'));
%     
% 
%     subplot(2,1,2);
%     plot(saveBinsFRS{1}(:,1),saveBinsFRS{1}(:,2))
%     hold on 
%     plot(saveBinsFRS{2}(:,1),saveBinsFRS{2}(:,2))
%     plot(saveBinsFRS{3}(:,1),saveBinsFRS{3}(:,2))
%      xlim([0 maxValPlot])
%      ylim([0 22])
%     legend(angles)
%     xlabel('Speed mm/sec')
%     ylabel('spikes/sec')

% Velocity

f(1) = figure(3); clf; 
    subplot(2,1,1);
    oppAngle_x = saveBinsVmVel{1}((saveBinsVmVel{1}(:,1) < maxValPlot & saveBinsVmVel{1}(:,1) > minValPlot),1);
    oppAngle_y = saveBinsVmVel{1}((saveBinsVmVel{1}(:,1) < maxValPlot & saveBinsVmVel{1}(:,1) > minValPlot),2);
    plot(oppAngle_x,oppAngle_y, 'color', [0.5 0.5 0.5])
    p(2,[1 2]) = polyfit(oppAngle_x,oppAngle_y, 1);  
    y1 = polyval(p(2,[1 2]), oppAngle_x);
    hold on 
    %plot(oppAngle_x,y1,'color', [0 0 0])
    
    prefAngle_x =  saveBinsVmVel{3}((saveBinsVmVel{3}(:,1) < maxValPlot & saveBinsVmVel{3}(:,1) > minValPlot),1);
    prefAngle_y = saveBinsVmVel{3}((saveBinsVmVel{3}(:,1) < maxValPlot & saveBinsVmVel{3}(:,1) > minValPlot),2);
    plot(prefAngle_x,prefAngle_y,'color', [1 0 0])
    p(1,[1 2]) = polyfit(prefAngle_x,prefAngle_y, 1);  
    y1 = polyval(p(1,[1 2]), prefAngle_x);
    %plot(prefAngle_x,y1,'color', [1 0 0])

    xlim([minValPlot+1 maxValPlot-1])
    xlabel('Velocity mm/sec')
    ylabel('Membrane Potential (mV)')
    set(gcf, 'color','w') 
    box off
    
    subplot(2,1,2);
    oppAngle_x = saveBinsFRVel{1}((saveBinsFRVel{1}(:,1) < maxValPlot & saveBinsFRVel{1}(:,1) > minValPlot),1);
    oppAngle_y = saveBinsFRVel{1}((saveBinsFRVel{1}(:,1) < maxValPlot & saveBinsFRVel{1}(:,1) > minValPlot),2);
    plot(oppAngle_x,oppAngle_y, 'color', [0.5 0.5 0.5])
    p(2,[3 4]) = polyfit(oppAngle_x,oppAngle_y, 1);  
    y1 = polyval(p(2,[3 4]), oppAngle_x);
    hold on 
   % plot(oppAngle_x,y1,'color', [0 0 0])
    
    prefAngle_x =  saveBinsFRVel{3}((saveBinsFRVel{3}(:,1) < maxValPlot & saveBinsFRVel{3}(:,1) > minValPlot),1);
    prefAngle_y = saveBinsFRVel{3}((saveBinsFRVel{3}(:,1) < maxValPlot & saveBinsFRVel{3}(:,1) > minValPlot),2);
    plot(prefAngle_x,prefAngle_y,'color', [1 0 0])
    p(1,[3 4]) = polyfit(prefAngle_x,prefAngle_y, 1);  
    y1 = polyval(p(1,[3 4]), prefAngle_x);
    %plot(prefAngle_x,y1,'color', [1 0 0])
 

    xlim([minValPlot+1 maxValPlot-1])
    xlabel('Velocity mm/sec')
    ylabel('Firing Rate (spikes/s)')
    set(gcf, 'color','w')  
    box off

keep = input('Save? ','s');

if strcmp(keep, 'y')
    cd(fileName) 
    savefig(f,'ActivityVsVelocity.fig')
    savefig(g,'ActivityVsSpeed.fig')
    cd('/Users/elenawesteinde/Documents/EphysCode') %'C:\Code\EphysCode'
end

end


