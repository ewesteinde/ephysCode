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

%% Make cosine basis function
% Bprs = struct;  % make new struct 'Bprs'
% Bprs.nBasis = 9;  % number of basis vectors
% Bprs.peakRange = [50, 550]; % location of 1st and last cosine peaks
% Bprs.timeRange = [-2000, 2000]; % 1st and last time bins for basis
% Bprs.dt = 1; % time bin size
% Bprs.logScaling = 'log';  % 'log' or 'linear'
% Bprs.logOffset = 1.25;
% alpha = 0.5; 
% 
% % currently log base w/ 5 bases returns best R^2 value but peak is near 0ms
% % lag. But resulting mean B*CosBasis looks the cleanest rather than
% % involving a lot of rapid oscillations and so better performance may be
% % due to less introduced noise, increasing number of bases didn't seem to
% % help much, may need to try a non symmetrical log basis set and/or how I'm
% % creating the symmetry as the narrow middle cosine curves aren't smooth or
% % increase alpha value towards lasso
% 
% % Make basis
% [cosBasis_log,tt] = makeRaisedCosBasis(Bprs);
% 
% if strcmp(Bprs.logScaling, 'log')
%     cosBasis_new = [];
%     cosBasis_shift = circshift(cosBasis_log, -50 ); %-Bprs.peakRange(1)
%     %neg_basis = circshift(flip(cosBasis_shift),-100);
%     neg_basis = flip(cosBasis_shift);
%     for time = 1:length(cosBasis_log(:,1))
%         if cosBasis_shift(time,1) > neg_basis(time,1)
%             cosBasis_new(time,1) = cosBasis_shift(time,1); 
%         else
%             cosBasis_new(time,1) = neg_basis(time,1);
%         end
%     end
% %     cosBasis(cosBasis(:,1) > 1) = 1; 
%     cosBasis_new = [cosBasis_new  cosBasis_shift(:,2:end)]; 
%     cosBasis_new = [cosBasis_new  neg_basis(:,2:end)];
% cosBasis_log = cosBasis_new(1001:3001,:);
% tt = tt(1001:3001);
% end
% 
% 
% Bprs = struct;  % make new struct 'Bprs'
% Bprs.nBasis = 9;  % number of basis vectors
% Bprs.peakRange = [-300, 300]; % location of 1st and last cosine peaks
% Bprs.timeRange = [-1000, 1000]; % 1st and last time bins for basis
% Bprs.dt = 1; % time bin size
% Bprs.logScaling = 'linear';  % 'log' or 'linear'
% 
% 
% % currently log base w/ 5 bases returns best R^2 value but peak is near 0ms
% % lag. But resulting mean B*CosBasis looks the cleanest rather than
% % involving a lot of rapid oscillations and so better performance may be
% % due to less introduced noise, increasing number of bases didn't seem to
% % help much, may need to try a non symmetrical log basis set and/or how I'm
% % creating the symmetry as the narrow middle cosine curves aren't smooth or
% % increase alpha value towards lasso
% 
% % Make basis
% [cosBasis_linear,tt] = makeRaisedCosBasis(Bprs);
% cosBasis = [];
% %cosBasis = [cosBasis_log cosBasis_linear];
% cosBasis = cosBasis_linear;
% 
% 
% % ---- make plots of basis and summed basis --------
% figure(20);clf; 
% subplot(211);
% plot(tt,cosBasis, 'linewidth', 2);
% set(gca,'tickdir', 'out');
% xlabel('time lag');  ylabel('basis function value');
% title('basis');
% 
% subplot(212);  % summed basis (to show tiling behavior: sums to constant)
% plot(tt,sum(cosBasis,2), 'linewidth', 2);
% set(gca,'tickdir', 'out');
% xlabel('time lag');  
% ylabel('sum');
% title('summed basis vectors');

%% make cosine basis function
Bprs = struct;  % make new struct 'Bprs'
Bprs.nBasis = 5;  % number of basis vectors
Bprs.peakRange = [50, 550]; % location of 1st and last cosine peaks
Bprs.timeRange = [-2000, 2000]; % 1st and last time bins for basis
Bprs.dt = 1; % time bin size
Bprs.logScaling = 'log';  % 'log' or 'linear'
Bprs.logOffset = 1.25;
alpha = 0.5; 

% currently log base w/ 5 bases returns best R^2 value but peak is near 0ms
% lag. But resulting mean B*CosBasis looks the cleanest rather than
% involving a lot of rapid oscillations and so better performance may be
% due to less introduced noise, increasing number of bases didn't seem to
% help much, may need to try a non symmetrical log basis set and/or how I'm
% creating the symmetry as the narrow middle cosine curves aren't smooth or
% increase alpha value towards lasso

% Make basis
[cosBasis,tt] = makeRaisedCosBasis(Bprs);

if strcmp(Bprs.logScaling, 'log')
    cosBasis_new = [];
    cosBasis_shift = circshift(cosBasis, -50 ); %-Bprs.peakRange(1)
    neg_basis = flip(cosBasis_shift); %circshift(flip(cosBasis_shift),-150);
    for time = 1:length(cosBasis(:,1))
        if cosBasis_shift(time,1) > neg_basis(time,1)
            cosBasis_new(time,1) = cosBasis_shift(time,1); 
        else
            cosBasis_new(time,1) = neg_basis(time,1);
        end
    end
%     cosBasis(cosBasis(:,1) > 1) = 1; 
    cosBasis_new = [cosBasis_new  cosBasis_shift(:,2:end)]; 
    cosBasis_new = [cosBasis_new  neg_basis(:,2:end)];
cosBasis = cosBasis_new(1001:3001,:);
tt = tt(1001:3001);
end


% cosBasis = resample(cosBasis,1,10);
% Bprs.dt = 1;
% [~,tt] = makeRaisedCosBasis(Bprs);

% ---- make plots of basis and summed basis --------
figure(20);clf; 
% subplot(211);
plot(tt,cosBasis, 'k','linewidth', 2);
set(gca,'tickdir', 'out');
xlabel('time lag');  ylabel('basis function value');
%title('basis');
box off 
set(gcf,'color',[1 1 1])
% 
% subplot(212);  % summed basis (to show tiling behavior: sums to constant)
% plot(tt,sum(cosBasis,2), 'linewidth', 2);
% set(gca,'tickdir', 'out');
% xlabel('time lag');  
% ylabel('sum');
% title('summed basis vectors');
%% fit linear model to each chunk
% number of features to fit coefficients to = number of cos basis functions

B_Vm = [];
B_fRate = [];
R2_fRate = []; 
R2_Vm = [];
count = 1; 

for chunk = 1:length(CL_chunks)
    
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


% make test chunk
%     X = [];
%     start = CL_chunks(1, 1);
%     finish = CL_chunks(1, 2);
%     
%     f = fieldnames(processed_behaviourData);
%     for k=1:numel(f)
%         chunk_behaviourData.(f{k}) = processed_behaviourData.(f{k})(start:finish); 
%     end
% 
% 
%     f = fieldnames(processed_trialData);
%     for k=1:numel(f)
%         chunk_trialData.(f{k}) = processed_trialData.(f{k})(start:finish); 
%     end
% transforming the original forward velocity trace w/ all the cosine basis functions, implementing lag & history
behaviour_data = chunk_behaviourData.vel_for;
X = [];
new_vf = [];
time_window = max(tt);
for basis = 1:size(cosBasis,2)
    new_vf = zeros(1,length(behaviour_data));
    for data = 1:length(behaviour_data)    
        if data < time_window + 1
            basis_idx = time_window + 1 - data + 1; 
            vf_basis = (behaviour_data(data)* cosBasis(basis_idx:end,basis))';
            new_vf(1: time_window + data) = new_vf(1: time_window + data) + vf_basis;
        elseif data > length(behaviour_data) - time_window
            basis_idx = data - (length(behaviour_data) - time_window); 
            vf_basis = (behaviour_data(data)* cosBasis(1:length(cosBasis) - basis_idx , basis))';
            new_vf(data - time_window:end) = new_vf(data - time_window:end) + vf_basis; 
        else 
            vf_basis = (behaviour_data(data)* cosBasis(:,basis))';
            new_vf(data - time_window:time_window + data) = new_vf(data - time_window:time_window + data) + vf_basis;
        end
    end
    X(:,basis) = new_vf; 
end


% debugging plot
% figure(22);clf;   
% 
%         yyaxis left
%         plot(chunk_behaviourData.time,new_vf, 'k') 
%         ylabel('Vm')
%         ylim([-(max(new_vf)) max(new_vf)])
%         %ylim([-70 -45])
%         hold on 
%         yyaxis right
%         plot(chunk_behaviourData.time, chunk_behaviourData.vel_for, 'r')
%         ylim([-(max(chunk_behaviourData.vel_for)) max(chunk_behaviourData.vel_for)])
%         ylabel('Vf mm/sec')
        
%

% feature scaling 
    for feature = 1:size(X,2)
        X(:,feature) = (X(:,feature) - min(X(:,feature)))/(max(X(:,feature)) - min(X(:,feature))); 
    end

% debugging plot
% figure(23);clf;   
% 
%         yyaxis left
%         plot(chunk_behaviourData.time,X(:,1), 'k') 
%         ylabel('Vm')
%         %ylim([-70 -45])
%         hold on 
%         yyaxis right
%         plot(chunk_behaviourData.time,X(:,9), 'r')
%         ylabel('Vf mm/sec')
%

%     y_fRate = chunk_trialData.fRate_sec(1000 : end - 1000)';
%     y_Vm = chunk_trialData.scaledOutput_down(1000 : end - 1000)'; 
    
    y_fRate = chunk_trialData.fRate_sec';
    y_Vm = chunk_trialData.scaledOutput_down'; 
   
    
        
    %X = X(1000:end-1000,:);
    
        [B_fRate_temp,FitInfo_fRate] = lasso(X, y_fRate,'NumLambda',5, 'Alpha', alpha) ;
        [B_Vm_temp,FitInfo_Vm] = lasso(X, y_Vm, 'NumLambda',5, 'Alpha', alpha) ;
        
        B_fRate(:,count) = B_fRate_temp(:,1); 
        B_Vm(:,count) = B_Vm_temp(:,1);
        
%         Intercept_fRate(:,chunk) = FitInfo_fRate.Intercept(1);
%         Intercept_Vm(:,chunk) = FitInfo_Vm.Intercept(1);

        yCalc_fRate = FitInfo_fRate.Intercept(1) + X*B_fRate(:,count);
        
%         +  B_fRate(1,chunk)*X(:,1) +  B_fRate(2,chunk)*X(:,2) +  B_fRate(3,chunk)*X(:,3) + ... 
%         B_fRate(4,chunk)*X(:,4) +  B_fRate(5,chunk)*X(:,5) +  B_fRate(6,chunk)*X(:,6) +  B_fRate(7,chunk)*X(:,7) +  ...
%         B_fRate(8,chunk)*X(:,8) +  B_fRate(9,chunk)*X(:,9);
%     
        yCalc_Vm = FitInfo_Vm.Intercept(1) + X*B_Vm(:,count);
        
%         B_Vm(1,chunk)*X(:,1) + B_Vm(2,chunk)*X(:,2) + B_Vm(3,chunk)*X(:,3) + ...
%         B_Vm(4,chunk)*X(:,4) + B_Vm(5,chunk)*X(:,5) + B_Vm(6,chunk)*X(:,6) + B_Vm(7,chunk)*X(:,7) + B_Vm(8,chunk)*X(:,8) + ...
%         B_Vm(9,chunk)*X(:,9); 

        R2_Vm_chunk =  1 - sum((y_Vm - yCalc_Vm).^2)/sum((y_Vm - mean(y_Vm)).^2);
        R2_fRate_chunk = 1 - sum((y_fRate - yCalc_fRate).^2)/sum((y_fRate - mean(y_fRate)).^2);
        
        R2_Vm(count) = R2_Vm_chunk;
        R2_fRate(count) = R2_fRate_chunk;
        
        count = count + 1; 
        
end   

figure(1);clf;
    subplot(2,1,1)
        plot(R2_Vm)
        text(count-1,max(R2_Vm),num2str(mean(R2_Vm)))
        box off
    subplot(2,1,2)
        plot(R2_fRate)
        text(count-1,max(R2_fRate),num2str(mean(R2_fRate)))
        box off
figure(2);clf;
    subplot(2,1,1)
    plot(tt,cosBasis*B_Vm,'Color',[0.75 0.75 0.75])
    hold on
    plot(tt,mean(cosBasis*B_Vm,2),'k')
    box off
    
subplot(2,1,2)
    plot(tt,cosBasis*B_fRate,'Color',[0.75 0.75 0.75])
    hold on
    plot(tt,mean(cosBasis*B_fRate,2),'k')
    box off
    
  
    
    % is R2 the best way to evaluate performance? Compare to squared error
    % loss & absolute error loss 
    
% yCalc_fRate = X*B_fRate + Intercept_fRate; 
% yCalc_Vm = X*B_Vm + Intercept_Vm;

% 
% for chunk = 1:size(B_Vm,2) 
%     
% %        Intercept_fRate(chunk) +  B_fRate(1,chunk)*X(:,1) +  B_fRate(2,chunk)*X(:,2) +  B_fRate(3,chunk)*X(:,3) + ... 
% %         B_fRate(4,chunk)*X(:,4) +  B_fRate(5,chunk)*X(:,5) +  B_fRate(6,chunk)*X(:,6) +  B_fRate(7,chunk)*X(:,7) +  ...
% %         B_fRate(8,chunk)*X(:,8) +  B_fRate(9,chunk)*X(:,9);
% % %     
% %         yCalc_Vm = Intercept_Vm(chunk) + B_Vm(1,chunk)*X(:,1) + B_Vm(2,chunk)*X(:,2) + B_Vm(3,chunk)*X(:,3) + ...
% %         B_Vm(4,chunk)*X(:,4) + B_Vm(5,chunk)*X(:,5) + B_Vm(6,chunk)*X(:,6) + B_Vm(7,chunk)*X(:,7) + B_Vm(8,chunk)*X(:,8) + ...
% %         B_Vm(9,chunk)*X(:,9); 
%     
%     yCalc_fRate = X*B_fRate(:,chunk)+ Intercept_fRate(chunk);
%     yCalc_Vm = X*B_Vm(:,chunk)+ Intercept_Vm(chunk);
%     R2_Vm(chunk) = 1 - sum((y_Vm - yCalc_Vm).^2)/sum((y_Vm - mean(y_Vm)).^2); 
%     R2_fRate(chunk) = 1 - sum((y_fRate - yCalc_fRate).^2)/sum((y_fRate - mean(y_fRate)).^2);
% end
%%

lags_neg = -200:50:0;
count = 1;
lag_labels_neg = {};
for labels = 1:length(lags_neg)
    lag_labels_neg{count} = num2str(lags_neg(labels));
    count = count + 1; 
end

lags_pos = 0:50:200;
count = 1;
lag_labels_pos = {};
for labels = 1:length(lags_pos)
    lag_labels_pos{count} = num2str(lags_pos(labels));
    count = count + 1; 
end



chunks = 1:length(B_fRate(1,:));
figure(55);clf;
subplot(2,2,1)
hold on
    for coeff = 1:length(X(1,:))
        plot(chunks,B_Vm(1:5,:))
        legend(lag_labels_neg)
    end    
subplot(2,2,3)
hold on
    for coeff = 1:length(X(1,:))
        plot(chunks,B_Vm(1:5,:))
        legend(lag_labels_neg)
    end
    
subplot(2,2,2)
hold on
    for coeff = 1:length(X(1,:))
        plot(chunks,R2_fRate)
        legend(lag_labels_neg)
    end
subplot(2,2,4)
hold on
    for coeff = 1:length(X(1,:))
        plot(chunks,R2_Vm)
        legend(lag_labels_neg)
    end
    

figure(56);clf;
subplot(2,2,1)
hold on
    for coeff = 1:length(X(1,:))
        plot(chunks,B_Vm(5:end,:))
        legend(lag_labels_pos)
    end    
subplot(2,2,3)
hold on
    for coeff = 1:length(X(1,:))
        plot(chunks,B_Vm(5:end,:))
        legend(lag_labels_pos)
    end
    
subplot(2,2,2)
hold on
    for coeff = 1:length(X(1,:))
        plot(chunks,R2_fRate)
        legend(lag_labels_pos)
    end
subplot(2,2,4)
hold on
    for coeff = 1:length(X(1,:))
        plot(chunks,R2_Vm)
        legend(lag_labels_pos)
    end
    
    %% sanity check, run clean simple code on one behavioural parameter at a time
    
    
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

    y_fRate = chunk_trialData.fRate_sec';
    y_Vm = chunk_trialData.scaledOutput_down'; 
   
    X = ((chunk_behaviourData.angle - min(chunk_behaviourData.angle))/(max(chunk_behaviourData.angle) - min(chunk_behaviourData.angle)))'; 
   
        [B_fRate_temp,FitInfo_fRate] = lasso(X, y_fRate,'NumLambda',5, 'Alpha', 0.5) ;
        [B_Vm_temp,FitInfo_Vm] = lasso(X, y_Vm, 'NumLambda',5, 'Alpha', 0.5) ;
        
        B_fRate(chunk) = B_fRate_temp(1); 
        B_Vm(chunk) = B_Vm_temp(1);
        
        yCalc_fRate = FitInfo_fRate.Intercept(1) +  B_fRate(chunk)*X; 
        yCalc_Vm = FitInfo_Vm.Intercept(1) + B_Vm(chunk)*X;  
        
        R2_Vm_chunk =  1 - sum((y_Vm - yCalc_Vm).^2)/sum((y_Vm - mean(y_Vm)).^2);
        R2_fRate_chunk = 1 - sum((y_fRate - yCalc_fRate).^2)/sum((y_fRate - mean(y_fRate)).^2);

        
        R2_Vm(chunk) = R2_Vm_chunk; 
        R2_fRate(chunk) = R2_fRate_chunk; 
       
end