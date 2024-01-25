
% %% load in experiment 
% 
% clear
% close all
% rootPath = '/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys'; % ; 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'%change dep on comp
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
% load('trialMeta.mat');
% load('trialData.mat');
% % openfig('CrossCorr.fig');
% % openfig('ActivityVsVelocity_trend.fig');
% 
% if isfield(trialMeta, 'notes')
%     disp(trialMeta.notes)
% end
% 
% cd('/Users/elenawesteinde/Documents/EphysCode/Analysis_code');


%%
% tStart = 1;
% tEnd = 1 ;
% prefHead = -80;
% range = 100; 
% headings = [wrapTo180(prefHead+180),wrapTo180(prefHead+90), prefHead];

%%
function [] = xcorr_headingBin(prefHead,range, processed_trialData, processed_behaviourData, fileName,iStart,iEnd)

headings = [wrapTo180(prefHead+180),wrapTo180(prefHead+90), prefHead];
t = 1; % run through each 5min trial

% Bin timepoints based on their headings

    fRate =  processed_trialData.fRate_sec(iStart:iEnd);
    Vm = processed_trialData.smooth_Vm(iStart:iEnd);
    angle = processed_behaviourData.angle(iStart:iEnd);
    vy = abs(processed_behaviourData.vel_yaw(iStart:iEnd));
    vf = processed_behaviourData.vel_for(iStart:iEnd);
    vs = abs(processed_behaviourData.vel_side(iStart:iEnd));
    speed = sqrt(vf.^2 + vs.^2);


    count = 0; 

    saveBinsVm = cell(1,3);
    saveBinsFR = cell(1,3); 
    saveBinsVf = cell(1,3);
    saveBinsVy = cell(1,3);
    saveBinsSp = cell(1,3);
    saveBinsIndex = cell(1,3);

    for head = headings
        VinHead = []; 
        VminHead = []; 
        fRateinHead = []; 
        index = [];
        count = count + 1;
        
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
            
            index = cat(2,index1, index2); 
        else

            lim1 = head-range/2;
            lim2 = head+range/2;

            if abs(lim1) > 180
                index = find(angle >= wrapTo180(lim1) | angle <= lim2);
                upperlim = lim2;
                lowerlim = wrapTo180(lim1);
            elseif abs(lim2) > 180
                index = find(angle >= lim1 | angle <= wrapTo180(lim2));
                upperlim = wrapTo180(lim2);
                lowerlim = lim1;
            else
                index = find(angle >= lim1 & angle <= lim2);
                upperlim = lim2;
                lowerlim = lim1;
            end
        end
        saveBinsIndex{count} = index;  
    end
    
 % extract continuous trial fragments within each heading bin 
     for b = 1:length(headings)
        count = 1;
        start = 0;
        frag = {};
        shiftPoint = [];
        for i = 2:length(saveBinsIndex{b})
            if saveBinsIndex{b}(i) - saveBinsIndex{b}(i-1) > 20
                shiftPoint(1,2) = i-1;
                if count == 1
                    frag{count} = saveBinsIndex{b}(1:shiftPoint(1,2)); 
                    count = count+1;
                    shiftPoint(1,1) = i;
                else
                    frag{count} = saveBinsIndex{b}(shiftPoint(1,1):shiftPoint(1,2));
                    count = count+1;
                    shiftPoint(1,1) = i;
                end 
            end

            if i == length(saveBinsIndex{b})
                if isempty(shiftPoint)
                    frag{count} = saveBinsIndex{b}(1:i);
                else
                    frag{count} = saveBinsIndex{b}(shiftPoint(1,1):i);
                end
            end
        end

        trial_fragments{b} = frag; 

% find trial fragments longer than one second & calculate the cross
% correlation b/w neural activity & velocity

% nanIDX  = find( isnan(velDot) );
% 
%     if ~isempty(nanIDX)
%         if nanIDX(1) == 1
%             velDot = velDot(nanIDX(end)+1:length(velDot)); 
%             Vm = Vm(nanIDX(end)+1:length(Vm)); 
%         elseif nanIDX(end) == length(velDot)
%             velDot = velDot(1:nanIDX(1)-1);
%             Vm = Vm(1:nanIDX(1)-1); 
%         else
%             error('NaN values inside velDot')
%         end
%     end
    
        count = 0; 
        xcorrSumVm_vf{b} = zeros(2001,1);
        xcorrSumFR_vf{b} = zeros(2001,1);
        
        xcorrSumVm_vy{b} = zeros(2001,1);
        xcorrSumFR_vy{b} = zeros(2001,1);
        
        xcorrSumVm_sp{b} = zeros(2001,1);
        xcorrSumFR_sp{b} = zeros(2001,1);
        
        
        for c = 1:length(trial_fragments{b})
            if length(trial_fragments{b}{c}) >= 4000

                vf_frag = vf(trial_fragments{b}{c}(1):trial_fragments{b}{c}(end));
                vy_frag = vy(trial_fragments{b}{c}(1):trial_fragments{b}{c}(end));
                sp_frag = speed(trial_fragments{b}{c}(1):trial_fragments{b}{c}(end));
                Vm_frag = Vm(trial_fragments{b}{c}(1):trial_fragments{b}{c}(end));
                FR_frag =  fRate(trial_fragments{b}{c}(1):trial_fragments{b}{c}(end));

                [cVm_vf, lagsVm_vf] = xcorr(Vm_frag',vf_frag,1000);
                [cFR_vf, lagsFR_vf] = xcorr(FR_frag',vf_frag,1000);
                [cVm_vy, lagsVm_vy] = xcorr(Vm_frag',vy_frag,1000);
                [cFR_vy, lagsFR_vy] = xcorr(FR_frag',vy_frag,1000);
                [cVm_sp, lagsVm_sp] = xcorr(Vm_frag',sp_frag,1000);
                [cFR_sp, lagsFR_sp] = xcorr(FR_frag',sp_frag,1000);

                % ,'normalized'
% 
%                 figure(1);clf
%                 subplot(2,1,1)
%                 plot(lagsVm_frag,cVm_frag)
%                 subplot(2,1,2)
%                 plot(lagsFR_frag,cFR_frag)
%                 
%                 figure(2); clf; 
%                 yyaxis left
%                 plot(velDot_frag)
%                 yyaxis right
%                 plot(Vm_frag)

                xcorrSumVm_vf{b} = xcorrSumVm_vf{b} + cVm_vf; 
                xcorrSumFR_vf{b} = xcorrSumFR_vf{b} + cFR_vf;
                
                xcorrSumVm_vy{b} = xcorrSumVm_vy{b} + cVm_vy; 
                xcorrSumFR_vy{b} = xcorrSumFR_vy{b} + cFR_vy;
                
                xcorrSumVm_sp{b} = xcorrSumVm_sp{b} + cVm_sp; 
                xcorrSumFR_sp{b} = xcorrSumFR_sp{b} + cFR_sp;

                count = count + 1;
            end
        end
         xcorrAveVm_vf{b}{1} = xcorrSumVm_vf{b}./count;
         xcorrAveVm_vf{b}{2} = lagsVm_vf;
         xcorrAveFR_vf{b}{1} = xcorrSumFR_vf{b}./count;
         xcorrAveFR_vf{b}{2} = lagsFR_vf;
         
         xcorrAveVm_vy{b}{1} = xcorrSumVm_vy{b}./count;
         xcorrAveVm_vy{b}{2} = lagsVm_vy;
         xcorrAveFR_vy{b}{1} = xcorrSumFR_vy{b}./count;
         xcorrAveFR_vy{b}{2} = lagsFR_vy;
         
         xcorrAveVm_sp{b}{1} = xcorrSumVm_sp{b}./count;
         xcorrAveVm_sp{b}{2} = lagsVm_sp;
         xcorrAveFR_sp{b}{1} = xcorrSumFR_sp{b}./count;
         xcorrAveFR_sp{b}{2} = lagsFR_sp;

%          figure(3);clf
%          subplot(2,1,1)
%          plot(xcorrAveVm{b}{2}, xcorrAveVm{b}{1})
%          if isempty(find( isnan(xcorrAveVm{b}{1})))
%             xline(xcorrAveVm{b}{2}(xcorrAveVm{b}{1} == max(xcorrAveVm{b}{1})),'-',num2str(xcorrAveVm{b}{2}(xcorrAveVm{b}{1} == max(xcorrAveVm{b}{1}))))
%          end
%          subplot(2,1,2)
%          plot(xcorrAveFR{b}{2}, xcorrAveFR{b}{1})
%          if isempty(find( isnan(xcorrAveVm{b}{1})))
%             xline(xcorrAveFR{b}{2}(xcorrAveFR{b}{1} == max(xcorrAveFR{b}{1})),'-',num2str(xcorrAveFR{b}{2}(xcorrAveFR{b}{1} == max(xcorrAveFR{b}{1}))))
%          end
         

     end % 1 heading 1 trial
    trialAve_headCorrVm_vf{t} = xcorrAveVm_vf;
    trialAve_headCorrFR_vf{t} = xcorrAveFR_vf;
    
    trialAve_headCorrVm_vy{t} = xcorrAveVm_vy;
    trialAve_headCorrFR_vy{t} = xcorrAveFR_vy;
    
    trialAve_headCorrVm_sp{t} = xcorrAveVm_sp;
    trialAve_headCorrFR_sp{t} = xcorrAveFR_sp;
% trialAve_headCorrVm{t}{b} t = trial num b = heading (1 = opp, 3 = pref)


% plot the data across trials & headings

% plot Vm corr
g = figure(); clf; 


subplot(3,3,1)
     plot(trialAve_headCorrVm_vf{1} {1}{2}, trialAve_headCorrVm_vf{1} {1}{1})
     if isempty(find( isnan(trialAve_headCorrVm_vf{1} {1}{1})))
        xline(trialAve_headCorrVm_vf{1} {1}{2}(trialAve_headCorrVm_vf{1} {1}{1} == max(trialAve_headCorrVm_vf{1} {1}{1})),'-',num2str(trialAve_headCorrVm_vf{1} {1}{2}(trialAve_headCorrVm_vf{1} {1}{1} == max(trialAve_headCorrVm_vf{1} {1}{1}))))
     end
    ylabel('For')
    title('opp head')
subplot(3,3,2)
     plot(trialAve_headCorrVm_vf{1} {2}{2}, trialAve_headCorrVm_vf{1} {2}{1})
     if isempty(find( isnan(trialAve_headCorrVm_vf{1} {2}{1})))
        xline(trialAve_headCorrVm_vf{1} {2}{2}(trialAve_headCorrVm_vf{1} {2}{1} == max(trialAve_headCorrVm_vf{1} {2}{1})),'-',num2str(trialAve_headCorrVm_vf{1} {2}{2}(trialAve_headCorrVm_vf{1} {2}{1} == max(trialAve_headCorrVm_vf{1} {2}{1}))))
     end
    title('int head')
subplot(3,3,3)
     plot(trialAve_headCorrVm_vf{1} {3}{2}, trialAve_headCorrVm_vf{1} {3}{1})
     if isempty(find( isnan(trialAve_headCorrVm_vf{1} {3}{1})))
        xline(trialAve_headCorrVm_vf{1} {3}{2}(trialAve_headCorrVm_vf{1} {3}{1} == max(trialAve_headCorrVm_vf{1} {3}{1})),'-',num2str(trialAve_headCorrVm_vf{1} {3}{2}(trialAve_headCorrVm_vf{1} {3}{1} == max(trialAve_headCorrVm_vf{1} {3}{1}))))
     end
        title('pref head')
    
subplot(3,3,4)
     plot(trialAve_headCorrVm_vy{1} {1}{2}, trialAve_headCorrVm_vy{1} {1}{1})
     if isempty(find( isnan(trialAve_headCorrVm_vy{1} {1}{1})))
        xline(trialAve_headCorrVm_vy{1} {1}{2}(trialAve_headCorrVm_vy{1} {1}{1} == max(trialAve_headCorrVm_vy{1} {1}{1})),'-',num2str(trialAve_headCorrVm_vy{1} {1}{2}(trialAve_headCorrVm_vy{1} {1}{1} == max(trialAve_headCorrVm_vy{1} {1}{1}))))
     end
     ylabel('Yaw')
subplot(3,3,5)
     plot(trialAve_headCorrVm_vy{1} {2}{2}, trialAve_headCorrVm_vy{1} {2}{1})
     if isempty(find( isnan(trialAve_headCorrVm_vy{1} {2}{1})))
        xline(trialAve_headCorrVm_vy{1} {2}{2}(trialAve_headCorrVm_vy{1} {2}{1} == max(trialAve_headCorrVm_vy{1} {2}{1})),'-',num2str(trialAve_headCorrVm_vy{1} {2}{2}(trialAve_headCorrVm_vy{1} {2}{1} == max(trialAve_headCorrVm_vy{1} {2}{1}))))
     end
subplot(3,3,6)
     plot(trialAve_headCorrVm_vy{1} {3}{2}, trialAve_headCorrVm_vy{1} {3}{1})
     if isempty(find( isnan(trialAve_headCorrVm_vy{1} {3}{1})))
        xline(trialAve_headCorrVm_vy{1} {3}{2}(trialAve_headCorrVm_vy{1} {3}{1} == max(trialAve_headCorrVm_vy{1} {3}{1})),'-',num2str(trialAve_headCorrVm_vy{1} {3}{2}(trialAve_headCorrVm_vy{1} {3}{1} == max(trialAve_headCorrVm_vy{1} {3}{1}))))
     end
     
 subplot(3,3,7)
     plot(trialAve_headCorrVm_sp{1} {1}{2}, trialAve_headCorrVm_sp{1} {1}{1})
     if isempty(find( isnan(trialAve_headCorrVm_sp{1} {1}{1})))
        xline(trialAve_headCorrVm_sp{1} {1}{2}(trialAve_headCorrVm_sp{1} {1}{1} == max(trialAve_headCorrVm_sp{1} {1}{1})),'-',num2str(trialAve_headCorrVm_sp{1} {1}{2}(trialAve_headCorrVm_sp{1} {1}{1} == max(trialAve_headCorrVm_sp{1} {1}{1}))))
     end
     ylabel('Speed')
subplot(3,3,8)
     plot(trialAve_headCorrVm_vy{1} {2}{2}, trialAve_headCorrVm_vy{1} {2}{1})
     if isempty(find( isnan(trialAve_headCorrVm_vy{1} {2}{1})))
        xline(trialAve_headCorrVm_sp{1} {2}{2}(trialAve_headCorrVm_sp{1} {2}{1} == max(trialAve_headCorrVm_sp{1} {2}{1})),'-',num2str(trialAve_headCorrVm_sp{1} {2}{2}(trialAve_headCorrVm_sp{1} {2}{1} == max(trialAve_headCorrVm_sp{1} {2}{1}))))
     end
subplot(3,3,9)
     plot(trialAve_headCorrVm_vy{1} {3}{2}, trialAve_headCorrVm_vy{1} {3}{1})
     if isempty(find( isnan(trialAve_headCorrVm_vy{1} {3}{1})))
        xline(trialAve_headCorrVm_sp{1} {3}{2}(trialAve_headCorrVm_sp{1} {3}{1} == max(trialAve_headCorrVm_sp{1} {3}{1})),'-',num2str(trialAve_headCorrVm_sp{1} {3}{2}(trialAve_headCorrVm_sp{1} {3}{1} == max(trialAve_headCorrVm_sp{1} {3}{1}))))
     end
     
    
% subplot(3,3,4)
%      plot(trialAve_headCorrVm{2} {1}{2}, trialAve_headCorrVm{2} {1}{1})
%      if isempty(find( isnan(trialAve_headCorrVm{2} {1}{1})))
%         xline(trialAve_headCorrVm{2} {1}{2}(trialAve_headCorrVm{2} {1}{1} == max(trialAve_headCorrVm{2} {1}{1})),'-',num2str(trialAve_headCorrVm{2} {1}{2}(trialAve_headCorrVm{2} {1}{1} == max(trialAve_headCorrVm{2} {1}{1}))))
%      end
%      ylabel('T1')
% subplot(3,3,5)
%      plot(trialAve_headCorrVm{2} {2}{2}, trialAve_headCorrVm{2} {2}{1})
%      if isempty(find( isnan(trialAve_headCorrVm{2} {2}{1})))
%         xline(trialAve_headCorrVm{2} {2}{2}(trialAve_headCorrVm{2} {2}{1} == max(trialAve_headCorrVm{2} {2}{1})),'-',num2str(trialAve_headCorrVm{2} {2}{2}(trialAve_headCorrVm{2} {2}{1} == max(trialAve_headCorrVm{2} {2}{1}))))
%      end
% subplot(3,3,6)
%      plot(trialAve_headCorrVm{2} {3}{2}, trialAve_headCorrVm{2} {3}{1})
%      if isempty(find( isnan(trialAve_headCorrVm{2} {3}{1})))
%         xline(trialAve_headCorrVm{2} {3}{2}(trialAve_headCorrVm{2} {3}{1} == max(trialAve_headCorrVm{2} {3}{1})),'-',num2str(trialAve_headCorrVm{2} {3}{2}(trialAve_headCorrVm{2} {3}{1} == max(trialAve_headCorrVm{2} {3}{1}))))
%      end
%      
%     subplot(3,3,7)
%      plot(trialAve_headCorrVm{3} {1}{2}, trialAve_headCorrVm{3} {1}{1})
%     if isempty(find( isnan(trialAve_headCorrVm{3} {1}{1})))
%         xline(trialAve_headCorrVm{3} {1}{2}(trialAve_headCorrVm{3} {1}{1} == max(trialAve_headCorrVm{3} {1}{1})),'-',num2str(trialAve_headCorrVm{3} {1}{2}(trialAve_headCorrVm{3} {1}{1} == max(trialAve_headCorrVm{3} {1}{1}))))
%     end
%      ylabel('T1')
% 
% subplot(3,3,8)
%      plot(trialAve_headCorrVm{3} {2}{2}, trialAve_headCorrVm{3} {2}{1})
%      if isempty(find( isnan(trialAve_headCorrVm{3} {2}{1})))
%         xline(trialAve_headCorrVm{3} {2}{2}(trialAve_headCorrVm{3} {2}{1} == max(trialAve_headCorrVm{3} {2}{1})),'-',num2str(trialAve_headCorrVm{3} {2}{2}(trialAve_headCorrVm{3} {2}{1} == max(trialAve_headCorrVm{3} {2}{1}))))
%      end
% 
% subplot(3,3,9)
%      plot(trialAve_headCorrVm{3} {3}{2}, trialAve_headCorrVm{3} {3}{1})
%      if isempty(find( isnan(trialAve_headCorrVm{3} {3}{1})))
%      xline(trialAve_headCorrVm{3} {3}{2}(trialAve_headCorrVm{3} {3}{1} == max(trialAve_headCorrVm{3} {3}{1})),'-',num2str(trialAve_headCorrVm{3} {3}{2}(trialAve_headCorrVm{3} {3}{1} == max(trialAve_headCorrVm{3} {3}{1}))))
%      end


    
    
 % plot FR corr across    
 
h = figure(); clf; 


subplot(3,3,1)
     plot(trialAve_headCorrFR_vf{1} {1}{2}, trialAve_headCorrFR_vf{1} {1}{1})
     if isempty(find( isnan(trialAve_headCorrFR_vf{1} {1}{1})))
        xline(trialAve_headCorrFR_vf{1} {1}{2}(trialAve_headCorrFR_vf{1} {1}{1} == max(trialAve_headCorrFR_vf{1} {1}{1})),'-',num2str(trialAve_headCorrFR_vf{1} {1}{2}(trialAve_headCorrFR_vf{1} {1}{1} == max(trialAve_headCorrFR_vf{1} {1}{1}))))
     end
    ylabel('For')
    title('opp head')
subplot(3,3,2)
     plot(trialAve_headCorrFR_vf{1} {2}{2}, trialAve_headCorrFR_vf{1} {2}{1})
     if isempty(find( isnan(trialAve_headCorrFR_vf{1} {2}{1})))
        xline(trialAve_headCorrFR_vf{1} {2}{2}(trialAve_headCorrFR_vf{1} {2}{1} == max(trialAve_headCorrFR_vf{1} {2}{1})),'-',num2str(trialAve_headCorrFR_vf{1} {2}{2}(trialAve_headCorrFR_vf{1} {2}{1} == max(trialAve_headCorrFR_vf{1} {2}{1}))))
     end
    title('int head')
subplot(3,3,3)
     plot(trialAve_headCorrFR_vf{1} {3}{2}, trialAve_headCorrFR_vf{1} {3}{1})
     if isempty(find( isnan(trialAve_headCorrFR_vf{1} {3}{1})))
        xline(trialAve_headCorrFR_vf{1} {3}{2}(trialAve_headCorrFR_vf{1} {3}{1} == max(trialAve_headCorrFR_vf{1} {3}{1})),'-',num2str(trialAve_headCorrFR_vf{1} {3}{2}(trialAve_headCorrFR_vf{1} {3}{1} == max(trialAve_headCorrFR_vf{1} {3}{1}))))
     end
        title('pref head')
    
subplot(3,3,4)
     plot(trialAve_headCorrFR_vy{1} {1}{2}, trialAve_headCorrFR_vy{1} {1}{1})
     if isempty(find( isnan(trialAve_headCorrFR_vy{1} {1}{1})))
        xline(trialAve_headCorrFR_vy{1} {1}{2}(trialAve_headCorrFR_vy{1} {1}{1} == max(trialAve_headCorrFR_vy{1} {1}{1})),'-',num2str(trialAve_headCorrFR_vy{1} {1}{2}(trialAve_headCorrFR_vy{1} {1}{1} == max(trialAve_headCorrFR_vy{1} {1}{1}))))
     end
     ylabel('Yaw')
subplot(3,3,5)
     plot(trialAve_headCorrFR_vy{1} {2}{2}, trialAve_headCorrFR_vy{1} {2}{1})
     if isempty(find( isnan(trialAve_headCorrFR_vy{1} {2}{1})))
        xline(trialAve_headCorrFR_vy{1} {2}{2}(trialAve_headCorrFR_vy{1} {2}{1} == max(trialAve_headCorrFR_vy{1} {2}{1})),'-',num2str(trialAve_headCorrFR_vy{1} {2}{2}(trialAve_headCorrFR_vy{1} {2}{1} == max(trialAve_headCorrFR_vy{1} {2}{1}))))
     end
subplot(3,3,6)
     plot(trialAve_headCorrFR_vy{1} {3}{2}, trialAve_headCorrFR_vy{1} {3}{1})
     if isempty(find( isnan(trialAve_headCorrFR_vy{1} {3}{1})))
        xline(trialAve_headCorrFR_vy{1} {3}{2}(trialAve_headCorrFR_vy{1} {3}{1} == max(trialAve_headCorrFR_vy{1} {3}{1})),'-',num2str(trialAve_headCorrFR_vy{1} {3}{2}(trialAve_headCorrFR_vy{1} {3}{1} == max(trialAve_headCorrFR_vy{1} {3}{1}))))
     end
     
subplot(3,3,7)
     plot(trialAve_headCorrFR_sp{1} {1}{2}, trialAve_headCorrFR_sp{1} {1}{1})
     if isempty(find( isnan(trialAve_headCorrFR_sp{1} {1}{1})))
        xline(trialAve_headCorrFR_sp{1} {1}{2}(trialAve_headCorrFR_sp{1} {1}{1} == max(trialAve_headCorrFR_sp{1} {1}{1})),'-',num2str(trialAve_headCorrFR_sp{1} {1}{2}(trialAve_headCorrFR_sp{1} {1}{1} == max(trialAve_headCorrFR_sp{1} {1}{1}))))
     end
     ylabel('Speed')
subplot(3,3,8)
     plot(trialAve_headCorrFR_sp{1} {2}{2}, trialAve_headCorrFR_sp{1} {2}{1})
     if isempty(find( isnan(trialAve_headCorrFR_sp{1} {2}{1})))
        xline(trialAve_headCorrFR_sp{1} {2}{2}(trialAve_headCorrFR_sp{1} {2}{1} == max(trialAve_headCorrFR_sp{1} {2}{1})),'-',num2str(trialAve_headCorrFR_sp{1} {2}{2}(trialAve_headCorrFR_sp{1} {2}{1} == max(trialAve_headCorrFR_sp{1} {2}{1}))))
     end
subplot(3,3,9)
     plot(trialAve_headCorrFR_sp{1} {3}{2}, trialAve_headCorrFR_sp{1} {3}{1})
     if isempty(find( isnan(trialAve_headCorrFR_sp{1} {3}{1})))
        xline(trialAve_headCorrFR_sp{1} {3}{2}(trialAve_headCorrFR_sp{1} {3}{1} == max(trialAve_headCorrFR_sp{1} {3}{1})),'-',num2str(trialAve_headCorrFR_sp{1} {3}{2}(trialAve_headCorrFR_sp{1} {3}{1} == max(trialAve_headCorrFR_sp{1} {3}{1}))))
     end
     
    
% subplot(3,3,4)
%      plot(trialAve_headCorrFR{2} {1}{2}, trialAve_headCorrFR{2} {1}{1})
%      if isempty(find( isnan(trialAve_headCorrFR{2} {1}{1})))
%         xline(trialAve_headCorrFR{2} {1}{2}(trialAve_headCorrFR{2} {1}{1} == max(trialAve_headCorrFR{2} {1}{1})),'-',num2str(trialAve_headCorrFR{2} {1}{2}(trialAve_headCorrFR{2} {1}{1} == max(trialAve_headCorrFR{2} {1}{1}))))
%      end
%      ylabel('T2')
% subplot(3,3,5)
%      plot(trialAve_headCorrFR{2} {2}{2}, trialAve_headCorrFR{2} {2}{1})
%      if isempty(find( isnan(trialAve_headCorrFR{2} {2}{1})))
%         xline(trialAve_headCorrFR{2} {2}{2}(trialAve_headCorrFR{2} {2}{1} == max(trialAve_headCorrFR{2} {2}{1})),'-',num2str(trialAve_headCorrFR{2} {2}{2}(trialAve_headCorrFR{2} {2}{1} == max(trialAve_headCorrFR{2} {2}{1}))))
%      end
% subplot(3,3,6)
%      plot(trialAve_headCorrFR{2} {3}{2}, trialAve_headCorrFR{2} {3}{1})
%      if isempty(find( isnan(trialAve_headCorrFR{2} {3}{1})))
%         xline(trialAve_headCorrFR{2} {3}{2}(trialAve_headCorrFR{2} {3}{1} == max(trialAve_headCorrFR{2} {3}{1})),'-',num2str(trialAve_headCorrFR{2} {3}{2}(trialAve_headCorrFR{2} {3}{1} == max(trialAve_headCorrFR{2} {3}{1}))))
%      end
%      
%     subplot(3,3,7)
%      plot(trialAve_headCorrFR{3} {1}{2}, trialAve_headCorrFR{3} {1}{1})
%     if isempty(find( isnan(trialAve_headCorrFR{3} {1}{1})))
%         xline(trialAve_headCorrFR{3} {1}{2}(trialAve_headCorrFR{3} {1}{1} == max(trialAve_headCorrFR{3} {1}{1})),'-',num2str(trialAve_headCorrFR{3} {1}{2}(trialAve_headCorrFR{3} {1}{1} == max(trialAve_headCorrFR{3} {1}{1}))))
%     end
%      ylabel('T3')
% subplot(3,3,8)
%      plot(trialAve_headCorrFR{3} {2}{2}, trialAve_headCorrFR{3} {2}{1})
%      if isempty(find( isnan(trialAve_headCorrFR{3} {2}{1})))
%         xline(trialAve_headCorrFR{3} {2}{2}(trialAve_headCorrFR{3} {2}{1} == max(trialAve_headCorrFR{3} {2}{1})),'-',num2str(trialAve_headCorrFR{3} {2}{2}(trialAve_headCorrFR{3} {2}{1} == max(trialAve_headCorrFR{3} {2}{1}))))
%      end
%   
% subplot(3,3,9)
%      plot(trialAve_headCorrFR{3} {3}{2}, trialAve_headCorrFR{3} {3}{1})
%      if isempty(find( isnan(trialAve_headCorrFR{3} {3}{1})))
%      xline(trialAve_headCorrFR{3} {3}{2}(trialAve_headCorrFR{3} {3}{1} == max(trialAve_headCorrFR{3} {3}{1})),'-',num2str(trialAve_headCorrFR{3} {3}{2}(trialAve_headCorrFR{3} {3}{1} == max(trialAve_headCorrFR{3} {3}{1}))))
%      end
  
keep = input('Save? ','s');

    if strcmp(keep, 'y')
        cd(fileName) 
        savefig(g,'CrossCorr_headingVm.fig')
        savefig(h,'CrossCorr_headingFR.fig')
        cd('/Users/elenawesteinde/Documents/EphysCode') %'/Users/elenawesteinde/Documents/EphysCode''C:\Code\EphysCode'
    end
end

    