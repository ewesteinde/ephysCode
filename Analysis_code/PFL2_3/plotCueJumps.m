function plotCueJumps(rootDir)

folders = get_folders_ephys(rootDir);

for f = 1:size(folders,1)
    %try
        folder = folders(f).folder; 
        if strcmp(folder(end),'.')
            folder = folder(1:end-2); 
        end

        processedDir = fullfile(folder,'processedData');
        load(fullfile(processedDir,'pro_behaviourData.mat'))
        load(fullfile(processedDir,'pro_trialData.mat'))

        numTrials = length(pro_behaviourData);

        %% visualize individual trials

        for t = 1:numTrials
            bData = pro_behaviourData{t};
            tData = pro_trialData{t};
            [rho, theta] = CalculateAverageHeading_ephys(bData,1.5, 'all');
            [jump_array, transition, yChannel] = detect_jumps_ephys(bData.frY, 10, 1000);
            
            total_mov_mm = abs(bData.vel_for) + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw))*4.5;
            no0vel_idx = find(total_mov_mm > 2);
            angle = bData.angle; 
            angle = angle(no0vel_idx); 
            
            try
                activity = tData.smooth_Vm;
            catch
                activity = tData.smoothVm;
            end
            activity = activity(no0vel_idx);
            
            edges_angle = [-180:30:180];
            [activityAngle, centers_angle, ~] = binData(activity, angle, edges_angle);
            
            prefHead_smooth = smoothdata(activityAngle,'gaussian',5);
            prefHead = centers_angle(prefHead_smooth == max(prefHead_smooth));
            
            for jump = 1:size(jump_array,1)
                jumpIdx = [jump_array(jump,1):jump_array(jump,3)];
                jumpTime = [-10:1/1000:10];
                if jumpIdx(end) > length(yChannel)
                    jumpIdx = jumpIdx(jumpIdx <= length(yChannel));
                    jumpTime = jumpTime(1:length(jumpIdx)); 
                end

                figure();
                set(gcf,'color','w')
                set(gcf,'renderer','painters')
                ax1 = subplot(4,1,1);
                plot(jumpTime,bData.angle(jumpIdx))
                ax2 = subplot(4,1,2);
                plot(jumpTime,tData.scaledOutput(jumpIdx))
                ax3 = subplot(4,1,3);
                plot(jumpTime, bData.vel_for(jumpIdx))
                ax4 = subplot(4,1,4);
                plot(jumpTime, bData.vel_yaw(jumpIdx))
                pause(0.5)
                linkaxes([ax1, ax2, ax3, ax4],'x')
                
                saveas(gcf, fullfile(folder,'figures',['jumpPlot_',num2str(jump),'.fig']))
            end
            
        end
   
end
end