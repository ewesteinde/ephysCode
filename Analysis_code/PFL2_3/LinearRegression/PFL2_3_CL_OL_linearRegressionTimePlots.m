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

cd('C:\Code\EphysCode'); 

if iscell(processed_behaviourData)
    processed_behaviourData = processed_behaviourData{1};
end

if iscell(processed_trialData)
    processed_trialData = processed_trialData{1};
end
%%

plotAngle_overlayVel_ScO(processed_trialData, processed_behaviourData)

%% Change if only part of the trial is usable
start =0;
finish = 2303;%processed_behaviourData.time(end);

iStart = find(processed_behaviourData.time == start); 
iEnd = find(processed_behaviourData.time == finish); 


f = fieldnames(processed_behaviourData);
for k=1:numel(f)
    processed_behaviourData.(f{k}) = processed_behaviourData.(f{k})(iStart:iEnd); 
end


f = fieldnames(processed_trialData);
for k=1:numel(f)
    processed_trialData.(f{k}) = processed_trialData.(f{k})(iStart:iEnd); 
end
       
%% separate idx at which fly is moving in CL vs OL
    
    startOL = 178.189*1000; % sec
    endOL = 203.249*1000;
    lengthCL = 179.99*1000;
    
% start = 1819.45;
% finish = 1844.47;
%     startOL  = find(processed_behaviourData.time == start); 
%     endOL = find(processed_behaviourData.time == finish); 

[CL, OL, CL_startStopIdx, OL_startStopIdx] = separateCL_OL(startOL, endOL, lengthCL, processed_behaviourData, processed_trialData);
%% CL_OL trial split into 1 min segements
% split 3 min continuous CL chunks into 1 min segments for more datapoints
count = 1;
CL_chunks = [];
for chunk = 1:length(CL_startStopIdx)
    distance = CL_startStopIdx(chunk,2) - CL_startStopIdx(chunk,1);
    step = floor(distance/3);
    l = length(CL_chunks); 
    CL_chunks(l+1,1) = CL_startStopIdx(chunk,1); 
    CL_chunks(l+1,2) = CL_startStopIdx(chunk,1) + (step);
    
    CL_chunks(l+2,1) = CL_startStopIdx(chunk,1) + (step)+1; 
    CL_chunks(l+2,2) = CL_startStopIdx(chunk,1) + (step)*2;
    
    CL_chunks(l+3,1) = CL_startStopIdx(chunk,1) + (step)*2+1; 
    CL_chunks(l+3,2) = CL_startStopIdx(chunk,2);
end

count = 1; 

%% Alternatively use a sliding 1 min window w/ overlap
CL_chunks = [];
window = 60; 
window = window * 1000;
step = 30 * 1000;
count = 1; 
for chunk = 1:length(CL_startStopIdx)
    chunk_idx = CL_startStopIdx(chunk, 1):CL_startStopIdx(chunk, 2);
    for i = 1:step:(length(chunk_idx) - window + 1)
        CL_chunks(count,1) = chunk_idx(i);
        CL_chunks(count,2) = chunk_idx(i + window - 1); 
        count = count + 1;
    end
    
    if (i + 2*window - step) > length(chunk_idx)
        CL_chunks(count,1) = chunk_idx(i + window - step);
        CL_chunks(count,2) = chunk_idx(end); 
        count = count + 1;
    end
end

%% jump/only CL trial split into 1 min segments

window = 60; 
window = window*1000; 
step = 30 * 1000; % determines degree of overlap if = to window then none 
CL_chunks = []; 
count = 1; 
for i = 1:step:length(processed_behaviourData.angle) - window + 1
    CL_chunks(count,1) = i; 
    CL_chunks(count,2) = i + window - 1; 
    count = count + 1; 
end

%%
% TO DO
    % repeat analysis with different lags b/w neural activity & heading/vel
    % likely need to implement cross validation 
    % look at combined parameters --> will need to then incorporate
    % regularization
    % develop smarter way to catch start/stop transition
    % develop more robust measurement of cell's pref head
    % test combined parameters (mult, add, etc)
    % figure out how to increase R2 for Frate
 
for chunk = 1:length(CL_chunks)
    X = [];
    start = CL_chunks(chunk, 1);
    finish = CL_chunks(chunk, 2);
    
    f = fieldnames(processed_behaviourData);
    for k=1:numel(f)
        chunk_behaviourData.(f{k}) = processed_behaviourData.(f{k})(start:finish); 
    end


    f = fieldnames(processed_trialData);
    for k=1:numel(f)
        chunk_trialData.(f{k}) = processed_trialData.(f{k})(start:finish); 
    end

    prefHeads = [-180:10:170];
    for p = 1:length(prefHeads)
        prefHead = prefHeads(p);
        
%     edges = [-180:30:180];
%     [centers_CL, mean_bin_CL] = create_binned_mean(chunk_behaviourData.angle, chunk_trialData.smooth_Vm, edges);
%     
%     prefHead(chunk) = centers_CL(mean_bin_CL == max(mean_bin_CL));
%     figure();
%     plot(centers_CL, mean_bin_CL,'-o');
%     xlabel('Angle');
%     ylabel('Vm'); 
%     title('CL heading pref')

% Need to take into account lag at some point soon
    y_Vm = chunk_trialData.scaledOutput_down';
        

        X(:,1) = (cosd(chunk_behaviourData.angle-prefHead) - min(cosd(chunk_behaviourData.angle-prefHead)))/(max(cosd(chunk_behaviourData.angle-prefHead))-min(cosd(chunk_behaviourData.angle-prefHead)));
%         X(:,1) = X_noLag((maxLag + 1 : end + minLag),1); 
        % implementing mean normalization
        X(:,2) = (chunk_behaviourData.vel_for - min(chunk_behaviourData.vel_for))/(max(chunk_behaviourData.vel_for) - min(chunk_behaviourData.vel_for)); 
%         X(:,2) = X_noLag((maxLag + 1 : end + minLag),2); 
        X(:,3) =(chunk_behaviourData.vel_yaw - min(chunk_behaviourData.vel_yaw))/(max(chunk_behaviourData.vel_yaw) - min(chunk_behaviourData.vel_yaw)); 

        
    % Calculate parameters

    %Matlab linear regression function
%         mdl_fRate = fitlm(X,y_fRate) ;
%         mdl_Vm = fitlm(X,y_Vm);
%         
%         
%         BTheta_fRate_chunk(p) = mdl_fRate.Coefficients{2,1};
%         Bvf_fRate_chunk(p) = mdl_fRate.Coefficients{3,1};
%         Bvy_fRate_chunk(p) = mdl_fRate.Coefficients{4,1};
% 
%         BTheta_Vm_chunk(p) = mdl_Vm.Coefficients{2,1}; 
%         Bvf_Vm_chunk(p) = mdl_Vm.Coefficients{3,1};
%         Bvy_Vm_chunk(p) = mdl_Vm.Coefficients{4,1};
%         
%         R2_Vm_chunk(p) = mdl_Vm.Rsquared.Ordinary;
%         R2_fRate_chunk(p) = mdl_fRate.Rsquared.Ordinary;
    % Lasso/elastic net regularization, must calc R^2 independently
        [B_Vm,FitInfo_Vm] = lasso(X, y_Vm, 'NumLambda',5, 'Alpha', 0.5) ;

        % Coefficients for heading & velocities

        BTheta_Vm_chunk(p) = B_Vm(1,1); 
        Bvf_Vm_chunk(p) = B_Vm(2,1);
        Bvy_Vm_chunk(p) = B_Vm(3,1);

        yCalc_Vm = FitInfo_Vm.Intercept(1) +  X*B_Vm(:,1); 
        
        R2_Vm_chunk(p) =  1 - sum((y_Vm - yCalc_Vm).^2)/sum((y_Vm - mean(y_Vm)).^2);
       

    end
    
    

    
    R2_Vm_posTheta = R2_Vm_chunk(BTheta_Vm_chunk >= 0);

    
    % set the chunk index to be the one where both the Bthetas are >= 0 and
    % R2 is the maximal value obtained across fitting the tested preferred headings. Controls for situations where two values of
    % preferred heading maximize R2 where they're 180 deg apart
    
    B_idx_Vm = BTheta_Vm_chunk >= 0 & R2_Vm_chunk == max(R2_Vm_posTheta); 
 
    
    R2_Vm(chunk) = R2_Vm_chunk(B_idx_Vm);
 

    BTheta_Vm(chunk) = BTheta_Vm_chunk(B_idx_Vm); 
    Bvf_Vm(chunk) = Bvf_Vm_chunk(B_idx_Vm);
    Bvy_Vm(chunk) = Bvy_Vm_chunk(B_idx_Vm);
    prefHead_Vm(chunk) = prefHeads(B_idx_Vm);

           
    [no0Vel_frag, no0Vel_idx] = remove0velocity(1.5, 1, chunk_behaviourData);
    %calculate percent time spent moving for each chunk
    percent_moving(chunk) = length(no0Vel_idx)/length(chunk_behaviourData.vel_for)*100;
    %calculate average euclidean speed for each chunk
    meanSpeed(chunk) = median(sqrt(chunk_behaviourData.vel_for(no0Vel_idx).^2 + chunk_behaviourData.vel_side(no0Vel_idx).^2));
    meanVf(chunk) = median(chunk_behaviourData.vel_for(no0Vel_idx)); 
    

%simply gather all timepoints during which the fly's speed was above the
%threshold

       f = fieldnames(chunk_behaviourData);
    for k=1:numel(f)
        no0vel_CL.(f{k}) = chunk_behaviourData.(f{k})(no0Vel_idx); 
    end
    
           f = fieldnames(chunk_trialData);
    for k=1:numel(f)
        no0vel_CL.(f{k}) = chunk_trialData.(f{k})(no0Vel_idx); 
    end
    
% calc heading vector & strength
    window = 60; %seconds
    minVel = 2; %mm/s
    sampRate = 1000; %Hz 

    [rho, theta, idx_windows, i, k] = PFL2_3_headingDist(window, minVel, chunk_behaviourData, sampRate, fileName, 0);
    
    % vector strength of heading
    ave_rho(chunk) = median(rho); 
    % circular mean of heading
    ave_theta(chunk) = median(theta); 
    ave_deg = rad2deg(ave_theta); 
    


end

%% run fit only on Vf data to find optimal lag

B_Vm = [];
B_fRate = []; 
for chunk = 1:length(CL_chunks)
    X = [];
    start = CL_chunks(chunk, 1);
    finish = CL_chunks(chunk, 2);
    
    f = fieldnames(processed_behaviourData);
    for k=1:numel(f)
        chunk_behaviourData.(f{k}) = processed_behaviourData.(f{k})(start:finish); 
    end


    f = fieldnames(processed_trialData);
    for k=1:numel(f)
        chunk_trialData.(f{k}) = processed_trialData.(f{k})(start:finish); 
    end
        

    minLag = -200;
    maxLag = 0;
    step = 50;

    y_fRate = chunk_trialData.fRate_sec(maxLag + 1 : end + minLag)';
    y_Vm = chunk_trialData.scaledOutput_down(maxLag + 1 : end + minLag)'; 
   
    X_noLag = (chunk_behaviourData.vel_for - min(chunk_behaviourData.vel_for))/(max(chunk_behaviourData.vel_for) - min(chunk_behaviourData.vel_for)); 
    count = 1;
    for lag = -200:50:0
        
        diffLag = minLag - lag; 
        mult = abs(diffLag/step); 
        X(:,1) = X_noLag((abs(minLag) - step * mult +1:end - step * mult));
        [B_fRate_temp,FitInfo_fRate] = lasso(X, y_fRate,'NumLambda',5, 'Alpha', 0.5) ;
        [B_Vm_temp,FitInfo_Vm] = lasso(X, y_Vm, 'NumLambda',5, 'Alpha', 0.5) ;
        
        B_fRate(count,chunk) = B_fRate_temp(1); 
        B_Vm(count,chunk) = B_Vm_temp(1);
        
        yCalc_fRate = FitInfo_fRate.Intercept(1) +  B_fRate(count,chunk)*X; 
        yCalc_Vm = FitInfo_Vm.Intercept(1) + B_Vm(count,chunk)*X;  
        
        R2_Vm_chunk =  1 - sum((y_Vm - yCalc_Vm).^2)/sum((y_Vm - mean(y_Vm)).^2);
        R2_fRate_chunk = 1 - sum((y_fRate - yCalc_fRate).^2)/sum((y_fRate - mean(y_fRate)).^2);

        
        R2_Vm(count,chunk) = R2_Vm_chunk; 
        R2_fRate(count,chunk) = R2_fRate_chunk; 
        
        count = count + 1;
        
    end            
end

%%
lags = -200:50:0;
count = 1;
lag_labels = {};
for labels = 1:length(lags)
    lag_labels{count} = num2str(lags(labels));
    count = count + 1; 
end
%%
chunks = 1:length(B_Vm(1,:));
figure(55);clf;
subplot(2,2,1)
hold on
    for coeff = 1:length(X(1,:))
        plot(chunks,B_Vm)
        legend(lag_labels)
    end    
subplot(2,2,3)
hold on
    for coeff = 1:length(X(1,:))
        plot(chunks,B_Vm)
        legend(lag_labels)
    end
    
subplot(2,2,2)
hold on
    for coeff = 1:length(X(1,:))
        plot(chunks,R2_fRate)
        legend(lag_labels)
    end
subplot(2,2,4)
hold on
    for coeff = 1:length(X(1,:))
        plot(chunks,R2_Vm)
        legend(lag_labels)
    end


%% Incorporating OL segments
if strcmp(trialMeta.trialType, 'CL_OL')
    for chunk = 1:length(OL_startStopIdx)
        X = [];
        start = OL_startStopIdx(chunk, 1);
        finish = OL_startStopIdx(chunk, 2);

        f = fieldnames(processed_behaviourData);
        for k=1:numel(f)
            chunk_behaviourData.(f{k}) = processed_behaviourData.(f{k})(start:finish); 
        end


        f = fieldnames(processed_trialData);
        for k=1:numel(f)
            chunk_trialData.(f{k}) = processed_trialData.(f{k})(start:finish); 
        end

        edges = [-180:30:180];
        [centers_CL, mean_bin_CL] = create_binned_mean(chunk_behaviourData.angle, chunk_trialData.smooth_Vm, edges);

        prefHead_OL(chunk) = centers_CL(mean_bin_CL == max(mean_bin_CL));
    %     figure();
    %     plot(centers_CL, mean_bin_CL,'-o');
    %     xlabel('Angle');
    %     ylabel('Vm'); 
    %     title('CL heading pref')

    % Need to take into account lag at some point soon

    %y = processed_trialData.smooth_Vm(start:finish)';
        y_fRate = chunk_trialData.fRate_sec';
        y_Vm = chunk_trialData.scaledOutput_down';

        X(:,1) = cosd(chunk_behaviourData.angle-prefHead_OL(chunk));
        X(:,2) = chunk_behaviourData.vel_for; 
        X(:,3) = chunk_behaviourData.vel_yaw; 
        %X(:,4) = sqrt(processed_behaviourData.vel_for(start:finish).^2 + processed_behaviourData.vel_side(start:finish).^2) ; 

        [m, n] = size(X);

    % Calculate parameters

    %Matlab linear regression function
        mdl_fRate = fitlm(X,y_fRate) ;
        mdl_Vm = fitlm(X,y_Vm);

        % Coefficients for heading & velocities
        BTheta_fRate_OL(chunk) = mdl_fRate.Coefficients{2,1};
        Bvf_fRate_OL(chunk) = mdl_fRate.Coefficients{3,1};
        Bvy_fRate_OL(chunk) = mdl_fRate.Coefficients{4,1};

        BTheta_Vm_OL(chunk) = mdl_Vm.Coefficients{2,1}; 
        Bvf_Vm_OL(chunk) = mdl_Vm.Coefficients{3,1};
        Bvy_Vm_OL(chunk) = mdl_Vm.Coefficients{4,1};


        % R^2
        R2_Vm_OL(chunk) = mdl_Vm.Rsquared.Ordinary;
        R2_fRate_OL(chunk) = mdl_fRate.Rsquared.Ordinary;


        [no0Vel_frag, no0Vel_idx] = remove0velocity(1.5, 1, chunk_behaviourData);
        %calculate percent time spent moving for each chunk
        percent_moving_OL(chunk) = length(no0Vel_idx)/length(chunk_behaviourData.vel_for)*100;
        %calculate average euclidean speed for each chunk
        meanSpeed_OL(chunk) = median(sqrt(chunk_behaviourData.vel_for(no0Vel_idx).^2 + chunk_behaviourData.vel_side(no0Vel_idx).^2));
        meanVf_OL(chunk) = median(chunk_behaviourData.vel_for(no0Vel_idx)); 

    end
OL_ID = [1:length(OL_startStopIdx)]*3; 
end
%%
%ave_deg
%ave_rho
%BTheta_fRate
%Bf_fRate
%By_fRate
%BTheta_Vm
%Bf_Vm
%By_Vm
%percent_moving
%meanSpeed
%prefHead

fig = figure(1);clf;
chunks = 1:length(BTheta_Vm);
subplot(7,1,1)
    plot(chunks,BTheta_Vm,'-or','MarkerFaceColor', 'r')
    ylabel('B Heading')
subplot(7,1,2);
    plot(chunks,Bvf_Vm,'-o','MarkerFaceColor',left_colour)
    ylabel('B Forward Vel')
%     hold on 
%     plot(OL_ID,Bvf_Vm_OL,'-ok','MarkerFaceColor', 'k')
subplot(7,1,3);
    plot(chunks,Bvy_Vm,'-og','MarkerFaceColor', 'g')
    ylabel('B Yaw Vel') 
subplot(7,1,4);
    plot(chunks,ave_deg,'-ok','MarkerFaceColor', 'k')
    ylabel('Heading')
    ylim([-200 200])
subplot(7,1,5);
    plot(chunks,ave_rho,'-ok','MarkerFaceColor', 'k')
    ylabel('Vector Strength')
    ylim([0 1])
subplot(7,1,6);
    plot(chunks,wrapTo180(prefHead_Vm),'-o','MarkerFaceColor',left_colour)
    ylim([-200 200])
    ylabel('Preferred Heading')
subplot(7,1,7)
    plot(chunks,R2_Vm,'-o','MarkerFaceColor',left_colour)
    ylabel('R^2 Vm')


%%
left_colour = [0 0.5 0];
right_colour = [0.5 0 0];

fig = figure(76);clf;
chunks = 1:length(BTheta_Vm);
name = strcat(date,', ',cell_num,', ',trial,' Linear Regression Summary Plot');
sgtitle(strcat(name));
g(1) = subplot(6,2,1);
    plot(chunks,BTheta_fRate,'-or','MarkerFaceColor', 'r')
    ylabel('BTh fRate')
%     hold on 
%     plot(OL_ID,BTheta_fRate_OL,'-ok','MarkerFaceColor', 'k')
g(2) = subplot(6,2,3);
    plot(chunks,Bvf_fRate,'-o','MarkerFaceColor',left_colour)
    ylabel('BVf fRate')
%     hold on 
%     plot(OL_ID,Bvf_fRate_OL,'-ok','MarkerFaceColor', 'k')
g(3) = subplot(6,2,5);
    plot(chunks,Bvy_fRate,'-og','MarkerFaceColor', 'g')
    ylabel('Bvy fRate')
%     hold on 
%     plot(OL_ID,Bvy_fRate_OL,'-ok','MarkerFaceColor', 'k')
g(4) = subplot(6,2,2);
    plot(chunks,BTheta_Vm,'-or','MarkerFaceColor', 'r')
    ylabel('BTh Vm')
%     hold on 
%     plot(OL_ID,BTheta_Vm_OL,'-ok','MarkerFaceColor', 'k')
g(5) = subplot(6,2,4);
    plot(chunks,Bvf_Vm,'-o','MarkerFaceColor',left_colour)
    ylabel('BVf Vm')
%     hold on 
%     plot(OL_ID,Bvf_Vm_OL,'-ok','MarkerFaceColor', 'k')
g(6) = subplot(6,2,6);
    plot(chunks,Bvy_Vm,'-og','MarkerFaceColor', 'g')
    ylabel('Bvy Vm')
%     hold on 
%     plot(OL_ID,Bvy_Vm_OL,'-ok','MarkerFaceColor', 'k')
g(7) = subplot(6,2,7);
    plot(chunks,ave_deg,'-ok','MarkerFaceColor', 'k')
    ylabel('heading')
    ylim([-200 200])
g(8) = subplot(6,2,8);
    plot(chunks,ave_rho,'-ok','MarkerFaceColor', 'k')
    ylabel('vector str')
    ylim([0 1])
g(9) = subplot(6,2,9);
    yyaxis left
        plot(chunks,wrapTo180(prefHead_Vm),'-o','MarkerFaceColor',left_colour)
    ylim([-200 200])
    ylabel('pHead Vm')
    yyaxis right
    plot(chunks,wrapTo180(prefHead_fRate),'-o','MarkerFaceColor',right_colour)
    ylim([-200 200])
    ylabel('pHead fRate')   
    set(fig,'defaultAxesColorOrder',[left_colour; right_colour]);    
%     hold on
%     plot(OL_ID, prefHead_OL,'b','MarkerFaceColor', 'b')
g(10) = subplot(6,2,10);
    miny = min(min(meanSpeed), min(meanVf));
    maxy = max(max(meanSpeed), max(meanVf));
    yyaxis left
    plot(chunks,meanSpeed,'-o','MarkerFaceColor', left_colour)
    ylabel('ave Speed')
    ylim([miny maxy])
    yyaxis right 
    plot(chunks,meanVf,'-o','MarkerFaceColor', right_colour)
    ylabel('vf')
    ylim([miny maxy])
g(11) = subplot(6,2,11);
    plot(chunks,percent_moving,'-o','MarkerFaceColor', 'k')
    ylabel('% moving')
    ylim([0 max(percent_moving)])
g(12) = subplot(6,2,12); 
    left_colour = [0 0.5 0];
    right_colour = [0.5 0 0];
    yyaxis left
    plot(chunks,R2_Vm,'-o','MarkerFaceColor',left_colour)
    ylabel('R^2 Vm')
    text(chunks(end)-1,max(R2_Vm),num2str(mean(R2_Vm)))
    yyaxis right
    plot(chunks,R2_fRate,'-o','MarkerFaceColor',right_colour)
    text(chunks(end)-1,max(R2_fRate),num2str(mean(R2_fRate)))
    ylabel('R^2 fRate')   
    set(fig,'defaultAxesColorOrder',[left_colour; right_colour]); 
linkaxes(g,'x')
%%
figure(55);clf;
            plot(chunks,Bvf_fRate_neg200)
            hold on
            plot(chunks,Bvf_fRate_neg150)
            plot(chunks,Bvf_fRate_neg100)
            plot(chunks,Bvf_fRate_neg50)
            plot(chunks,Bvf_fRate_pos50)
            plot(chunks,Bvf_fRate_pos100)
            plot(chunks,Bvf_fRate_pos150)
            plot(chunks,Bvf_fRate_pos200)
            plot(chunks,Bvf_fRate,'-ob','MarkerFaceColor', 'b') 
            legend 