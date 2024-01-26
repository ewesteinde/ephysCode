function [xPos, yPos] = plot2DTrajectory_activity(rootDir, savePlot)
    
folders = get_folders_ephys_behaviour(rootDir, 1); 

%% Process each folder
folderNum = length(folders);
fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
countFail = 1;


for ff = 1:folderNum
      %try
      %% Get folder information
        folder = folders(ff).folder;
        disp(folder) 
        load(fullfile(folder,'pro_behaviourData.mat'))
        load(fullfile(folder,'pro_trialData.mat'))
        
        bData = processed_behaviourData{1}; 
        tData = processed_trialData{1}; 
        
        newSampRate = 60; 

        vq = bData.time(1):1/newSampRate:bData.time(end);
        fn = fieldnames(bData);
        if istable(bData)
            fn = fn(1:end - 3); 
        end
        for field = 1:numel(fn)
            if isempty(regexp(fn{field},'disp'))
                bData_interp.(fn{field}) = interp1(bData.time,bData.(fn{field}),vq,'linear');
            end
        end  
        
        
        vq = tData.time(1):1/newSampRate:tData.time(end);
        fn = fieldnames(tData);
        if istable(tData)
            fn = fn(1:end - 3); 
        end
        for field = 1:numel(fn)
            if isempty(regexp(fn{field},'disp'))
                tData_interp.(fn{field}) = interp1(tData.time,tData.(fn{field}),vq,'linear');
            end
        end  
            %%

            %% Helen's path code

            %for nTrial = numTrials
                notNan_idx = find(~isnan(bData_interp.angle) & ~isnan(bData_interp.vel_for) & ~isnan(bData_interp.vel_side));
                yawAngPos = bData_interp.angle(notNan_idx);
                fwdAngVel = bData_interp.vel_for(notNan_idx);
                slideAngVel = bData_interp.vel_side(notNan_idx);
                time = bData_interp.time(notNan_idx);
                
                fRate = tData_interp.fRate_sec(notNan_idx);
                VM = tData_interp.smoothVm(notNan_idx);

                sampRate = newSampRate; 
                % conversion factor between degrees and mm
                circum = 9 * pi; % circumference of ball, in mm
                mmPerDeg = circum / 360; % mm per degree of ball

                % position incorporating heading - as if fly were walking on x-y plane,
                %  x-y coordinates at each time point
                % start with fly at (0,0) and facing 0 deg
                zeroedYawAngPos = yawAngPos - yawAngPos(1); 

                % movement in x (in degrees) at each time point
                xChangePos = (fwdAngVel ./ sampRate) .* sind(zeroedYawAngPos) + ...
                    (slideAngVel ./ sampRate) .* sind(zeroedYawAngPos + 90);  

                % x position in mm (i.e. x-coordinate of fly's position at each time 
                %  point), starts at 0
                xPos = (cumsum(xChangePos) - xChangePos(1)) .* mmPerDeg;
                minX = min(xPos);
                maxX = max(xPos);

                % movement in y (in degrees) at each time point
                yChangePos = (fwdAngVel ./ sampRate) .* cosd(zeroedYawAngPos) + ...
                    (slideAngVel ./ sampRate) .* cosd(zeroedYawAngPos + 90);

                % y position in mm (i.e. y-coordinate of fly's position at each time 
                %  point), starts at 0
                yPos = (cumsum(yChangePos) - yChangePos(1)) .* mmPerDeg;
                minY = min(yPos);
                maxY = max(yPos);
                
                VM_smooth = smoothdata(VM,'movmedian', 20) ;
                fwdAngVel_smooth = smoothdata(fwdAngVel,'movmedian', 20) ;

                figure();patch(xPos,yPos,VM_smooth,'EdgeColor','interp','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',5);
                %figure();patch('xData',xPos,'yData',yPos,'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
                hold on
                %patch(xPos(transition == 1),yPos(transition == 1),time(transition == 1),'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',7);
                %patch('xData',xPos(transition == 1),'yData',yPos(transition == 1),'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',7);
                colormap(jet)
                a=colorbar; a.Label.String = 'Vm';
                set(gcf,'Renderer','painters')
                set(gcf,'color','w');
                xlabel('mm')
                xlim([min(minX,minY),max(maxX, maxY)])
                ylim([min(minX,minY),max(maxX, maxY)])

                idxes = regexp(folder,'\');
                title(folder(idxes(end) + 1:end), 'Interpreter', 'none')
                
                if savePlot
                    saveas(gcf, fullfile(folder,'figures','pathTrajectory_Vm.fig'));
                end
                
                FR_smooth = smoothdata(fRate,'movmedian', 20) ;

                figure();patch(xPos,yPos,FR_smooth,'EdgeColor','interp','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',5);
                %figure();patch('xData',xPos,'yData',yPos,'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
                hold on
                %patch(xPos(transition == 1),yPos(transition == 1),time(transition == 1),'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',7);
                %patch('xData',xPos(transition == 1),'yData',yPos(transition == 1),'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',7);
                colormap(jet)
                a=colorbar; a.Label.String = 'FR';
                set(gcf,'Renderer','painters')
                set(gcf,'color','w');
                xlabel('mm')
                xlim([min(minX,minY),max(maxX, maxY)])
                ylim([min(minX,minY),max(maxX, maxY)])

                idxes = regexp(folder,'\');
                title(folder(idxes(end) + 1:end), 'Interpreter', 'none')
                
                if savePlot
                    saveas(gcf, fullfile(folder,'figures','pathTrajectory_FR.fig'));
                end
                
                figure();patch(xPos,yPos,fwdAngVel_smooth,'EdgeColor','interp','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',5);
                %figure();patch('xData',xPos,'yData',yPos,'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
                hold on
                %patch(xPos(transition == 1),yPos(transition == 1),time(transition == 1),'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',7);
                %patch('xData',xPos(transition == 1),'yData',yPos(transition == 1),'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',7);
                colormap(jet)
                a=colorbar; a.Label.String = 'vf';
                set(gcf,'Renderer','painters')
                set(gcf,'color','w');
                xlabel('mm')
                xlim([min(minX,minY),max(maxX, maxY)])
                ylim([min(minX,minY),max(maxX, maxY)])

                idxes = regexp(folder,'\');
                title(folder(idxes(end) + 1:end), 'Interpreter', 'none')
                
                if savePlot
                    saveas(gcf, fullfile(folder,'figures','pathTrajectory_vf.fig'));
                end


%                 figure();patch(ftData_dat.intX{1},ftData_dat.intY{1},ftData_dat.trialTime{1},'EdgeColor','interp','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
%                 a=colorbar; a.Label.String = 'Time (sec)';
%                 set(gcf,'color','w');
%                 axis off

%                 saveas(gcf, fullfile(folder,'images',[expID,'_',num2str(nTrial),'_pathTrajectory.fig']));
%                 saveas(gcf, fullfile(folder,'images',[expID,'_',num2str(nTrial),'_pathTrajectory.png']));
%                 saveas(gcf, fullfile(folder,'images',[expID,'_',num2str(nTrial),'_pathTrajectory.svg']));
            %end
%       catch
%           disp(['Folder: ',folder,' failed'])
%       end
end
end
%close all 