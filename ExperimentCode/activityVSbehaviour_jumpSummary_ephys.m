function activityVSbehaviour_jumpSummary_ephys(rootDir)       

    folders = get_folders_ephys(rootDir);
    countsum90 = 1; 
    countsum180 = 1; 
    G90count = 1; 
    G180count = 1;

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
        %try
        bData = pro_behaviourData{t};
        tData = pro_trialData{t};
        [rho, theta] = CalculateAverageHeading_ephys(bData,1.5, 'all');
        [jump_array, transition, yChannel] = detect_jumps_ephys(bData.frY, 10, 1000);

        total_mov_mm = abs(bData.vel_for) + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw))*4.5;
        no0vel_idx = find(total_mov_mm > 2);
        angle = bData.angle;

                count90 = 1;
                count180 = 1;
                jump90_time = []; 
                jump180_time = []; 
                jumps_90 = []; 
                jumps_180 = []; 
                for jump = 1:size(jump_array,1)
                    if abs(jump_array(jump,4)) == 90
                        jumps_90(count90,:) = [jump_array(jump,1):jump_array(jump,3)];
                        jump90_time(:,count90) = (jumps_90(count90,:) - jump_array(jump,2))/1000;
                        count90 = count90 + 1; 
                    else
                        jumps_180(count180,:) = [jump_array(jump,1):jump_array(jump,3)];
                        jump180_time(:,count180) = (jumps_180(count180,:) - jump_array(jump,2))/1000;
                        count180 = count180 + 1; 
                    end
                end
                
                lastIdx_preJump = (length(jump180_time) - 1)/2;
                for jump = 1:size(jumps_90,1)
                    try
                        diff_goal = rad2deg(angdiff(deg2rad(angle(jumps_90(jump,lastIdx_preJump))),theta)); 
                    catch
                        disp('what')
                    end

                    vf90(countsum90,:) = bData.vel_for(jumps_90(jump,:)); 
                    vy90(countsum90,:) = abs(bData.vel_yaw(jumps_90(jump,:))); 
                    vs90(countsum90,:) = abs(bData.vel_side(jumps_90(jump,:))); 
                    angle90(countsum90,:) = bData.angle(jumps_90(jump,:)); 
                    countsum90 = countsum90 + 1; 

                    if abs(diff_goal) < 30
                        jump90_meno(G90count).folder = folder;
                        jump90_meno(G90count).goal = theta; 
                        jump90_meno(G90count).rho = rho; 
                        jump90_meno(G90count).jumpIdx = jumps_90(jump,:);
                        G90count = G90count + 1; 
                    end

                end

                for jump = 1:size(jumps_180,1)
                    try
                        diff_goal = rad2deg(angdiff(deg2rad(angle(jumps_180(jump,lastIdx_preJump))),theta)); 
                
                        vf180(countsum180,:) = bData.vel_for(jumps_180(jump,:)); 
                        vy180(countsum180,:) = abs(bData.vel_yaw(jumps_180(jump,:))); 
                        vs180(countsum180,:) = abs(bData.vel_side(jumps_180(jump,:))); 
                        angle180(countsum180,:) = bData.angle(jumps_180(jump,:)); 
                        countsum180 = countsum180 + 1; 
                    catch
                        disp('idx exceeds count')
                    end

                    if abs(diff_goal) < 30
                        jump180_meno(G90count).folder = folder;
                        jump180_meno(G90count).goal = theta; 
                        jump180_meno(G90count).rho = rho; 
                        jump180_meno(G90count).jumpIdx = jumps_90(jump,:);
                        G180count = G180count + 1; 
                    end
                end

            end

%         catch
%             disp([folder, ' failed'])
%         end
    end  
    jump90_time = repmat(jump90_time(:,1),1,size(vf90,1)); 
    jump180_time = repmat(jump180_time(:,1),1,size(vf180,1));
    figure(); 
            ax1= subplot(3,1,1);
            set(gcf, 'color','w')
            set(gcf,'renderer','painters')
            plot(jump90_time',vf90);
            hold on
            plot(jump90_time(:,1)', mean(vf90,1,'omitnan'),'k','lineWidth',1.5)
            xline(0)
            ylabel('vf')
            box off

            ax2= subplot(3,1,2);
            plot(jump90_time,rad2deg(vy90)');
            hold on
            plot(jump90_time(:,1), mean(rad2deg(vy90),1,'omitnan'),'k','lineWidth',1.5)
            xline(0)
            ylabel('vy')
            box off

            ax3= subplot(3,1,3);
            plot(jump90_time,vs90');
            hold on
            plot(jump90_time(:,1), mean(vs90,1,'omitnan'),'k','lineWidth',1.5)
            xline(0)
            ylabel('vs')
            box off



            figure(); 
            ax1= subplot(3,1,1);
            set(gcf, 'color','w')
            set(gcf,'renderer','painters')
            plot(jump180_time,vf180');
            hold on
            plot(jump180_time(:,1), mean(vf180,1,'omitnan'),'k','lineWidth',1.5)
            xline(0)
            ylabel('vf')
            box off

            ax2= subplot(3,1,2);
            plot(jump180_time,rad2deg(vy180'));
            hold on
            plot(jump180_time(:,1), mean(rad2deg(vy180),1,'omitnan'),'k','lineWidth',1.5)
            xline(0)
            ylabel('vy')
            box off

            ax3= subplot(3,1,3);
            plot(jump180_time,vs180');
            hold on
            plot(jump180_time(:,1), mean(vs180,1,'omitnan'),'k','lineWidth',1.5)
            xline(0)
            box off

end