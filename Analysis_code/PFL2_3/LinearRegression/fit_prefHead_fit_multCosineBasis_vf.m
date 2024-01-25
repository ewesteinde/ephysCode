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
plot(tt,cosBasis, 'k','linewidth', 1);
set(gca,'tickdir', 'out');
xlabel('time lag');  ylabel('basis function value');
title('basis');
box off
set(gcf,'color',[1 1 1]); 

subplot(212);  % summed basis (to show tiling behavior: sums to constant)
plot(tt,sum(cosBasis,2), 'linewidth', 2);
set(gca,'tickdir', 'out');
xlabel('time lag');  
ylabel('sum');
title('summed basis vectors');


%%

tic

prefHeads = [-180:10:170];
B_Vm = [];
B_fRate = [];
R2_fRate = []; 
R2_Vm = [];

BTheta_Vm_chunk = zeros(1,length(prefHeads));
R2_Vm_chunk = zeros(1,length(prefHeads));
R2_Vm_heading = zeros(1,length(CL_chunks));
BTheta_Vm_heading = zeros(1,length(CL_chunks)); 
prefHead_Vm = zeros(1,length(length(CL_chunks)));

count = 1; 

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

y_fRate = log(chunk_trialData.fRate_sec');
y_fRate(y_fRate == -Inf) = 0;
y_Vm = chunk_trialData.scaledOutput_down';



for p = 1:length(prefHeads)
    
    prefHead = prefHeads(p);
    
    X = (cosd(chunk_behaviourData.angle-prefHead) - min(cosd(chunk_behaviourData.angle-prefHead)))/(max(cosd(chunk_behaviourData.angle-prefHead))-min(cosd(chunk_behaviourData.angle-prefHead)));
    X = X';
    
%     new_vf = [];
%     time_window = max(tt);
%     X_zero = zeros(length(chunk_behaviourData.vel_for), size(cosBasis,2));
%     X = [X X_zero];
%     for basis = 1:size(cosBasis,2)
%         new_vf = zeros(1,length(chunk_behaviourData.vel_for));
%         for data = 1:length(chunk_behaviourData.vel_for)    
%             if data < time_window + 1
%                 basis_idx = time_window + 1 - data + 1; 
%                 vf_basis = (chunk_behaviourData.vel_for(data)* cosBasis(basis_idx:end,basis))';
%                 new_vf(1: time_window + data) = new_vf(1: time_window + data) + vf_basis;
%             elseif data > length(chunk_behaviourData.vel_for) - time_window
%                 basis_idx = data - (length(chunk_behaviourData.vel_for) - time_window); 
%                 vf_basis = (chunk_behaviourData.vel_for(data)* cosBasis(1:length(cosBasis) - basis_idx , basis))';
%                 new_vf(data - time_window:end) = new_vf(data - time_window:end) + vf_basis; 
%             else 
%                 vf_basis = (chunk_behaviourData.vel_for(data)* cosBasis(:,basis))';
%                 new_vf(data - time_window:time_window + data) = new_vf(data - time_window:time_window + data) + vf_basis;
%             end
%         end
%         % feature scaling 
%         new_vf = (new_vf - min(new_vf))/(max(new_vf) - min(new_vf)); 
%         X(:,basis + 1) = new_vf; 
%     end
    
    [B_Vm,FitInfo_Vm] = lasso(X, y_Vm, 'NumLambda',2, 'Alpha', 0.5) ; 
    B_Vm_chunk(:,p) = B_Vm(:,1); 
    
    [B_fRate,FitInfo_fRate] = lasso(X, y_fRate, 'NumLambda',2, 'Alpha', 0.5) ; 
    B_fRate_chunk(:,p) = B_fRate(:,1);

    yCalc_Vm = FitInfo_Vm.Intercept(1) +  X*B_Vm(:,1);
    yCalc_fRate = FitInfo_fRate.Intercept(1) +  X*B_fRate(:,1);

    R2_Vm_chunk(p) =  1 - sum((y_Vm - yCalc_Vm).^2)/sum((y_Vm - mean(y_Vm)).^2);
    R2_fRate_chunk(p) =  1 - sum((y_fRate - yCalc_fRate).^2)/sum((y_fRate - mean(y_fRate)).^2);
end
    
    R2_Vm_posTheta = R2_Vm_chunk(B_Vm_chunk(1,:) >= 0);
    R2_fRate_posTheta = R2_fRate_chunk(B_fRate_chunk(1,:) >= 0);
    
    % set the chunk index to be the one where both the Bthetas are >= 0 and
    % R2 is the maximal value obtained across fitting the tested preferred headings. Controls for situations where two values of
    % preferred heading maximize R2 where they're 180 deg apart
    
    B_idx_Vm = B_Vm_chunk(1,:) >= 0 & R2_Vm_chunk == max(R2_Vm_posTheta); 
    B_idx_fRate = B_fRate_chunk(1,:) >= 0 & R2_fRate_chunk == max(R2_fRate_posTheta);
    
    R2_Vm_heading(chunk) = R2_Vm_chunk(B_idx_Vm);
    R2_fRate_heading(chunk) = R2_fRate_chunk(B_idx_fRate);

%     B_Vm(:,chunk) = B_Vm_chunk(:,B_idx_Vm); 
%     B_fRate(:,chunk) = B_fRate_chunk(:,B_idx_Vm);
    
    prefHead_Vm(chunk) = prefHeads(B_idx_Vm);
    prefHead_fRate(chunk) = prefHeads(B_idx_fRate);

% transforming the original forward velocity trace w/ all the cosine basis functions, implementing lag & history
% 
    new_vf = [];
    time_window = max(tt);
    X_fRate = [];
    X_Vm = [];
    X_zero = zeros(length(chunk_behaviourData.vel_for), size(cosBasis,2));
    X_Vm(:,1) = (cosd(chunk_behaviourData.angle-prefHead_Vm(chunk)) - min(cosd(chunk_behaviourData.angle-prefHead_Vm(chunk))))/(max(cosd(chunk_behaviourData.angle-prefHead_Vm(chunk)))-min(cosd(chunk_behaviourData.angle-prefHead_Vm(chunk))));
    X_fRate(:,1) = (cosd(chunk_behaviourData.angle-prefHead_fRate(chunk)) - min(cosd(chunk_behaviourData.angle-prefHead_fRate(chunk))))/(max(cosd(chunk_behaviourData.angle-prefHead_fRate(chunk)))-min(cosd(chunk_behaviourData.angle-prefHead_fRate(chunk))));
    X_Vm = [X_Vm X_zero];
    X_fRate = [X_fRate X_zero];
    
    for basis = 1:size(cosBasis,2)
        new_vf = zeros(1,length(chunk_behaviourData.vel_for));
        for data = 1:length(chunk_behaviourData.vel_for)    
            if data < time_window + 1
                basis_idx = time_window + 1 - data + 1; 
                vf_basis = (chunk_behaviourData.vel_for(data)* cosBasis(basis_idx:end,basis))';
                new_vf(1: time_window + data) = new_vf(1: time_window + data) + vf_basis;
            elseif data > length(chunk_behaviourData.vel_for) - time_window
                basis_idx = data - (length(chunk_behaviourData.vel_for) - time_window); 
                vf_basis = (chunk_behaviourData.vel_for(data)* cosBasis(1:length(cosBasis) - basis_idx , basis))';
                new_vf(data - time_window:end) = new_vf(data - time_window:end) + vf_basis; 
            else 
                vf_basis = (chunk_behaviourData.vel_for(data)* cosBasis(:,basis))';
                new_vf(data - time_window:time_window + data) = new_vf(data - time_window:time_window + data) + vf_basis;
            end
        end
        % feature scaling 
        new_vf = (new_vf - min(new_vf))/(max(new_vf) - min(new_vf)); 
        X_Vm(:,basis + 1) = new_vf; 
        X_fRate(:,basis + 1) = new_vf;
    end
    
%     X(:,2) = (chunk_behaviourData.vel_for - min(chunk_behaviourData.vel_for))/(max(chunk_behaviourData.vel_for) - min(chunk_behaviourData.vel_for));
%     X(:,3) = (chunk_behaviourData.vel_for - min(chunk_behaviourData.vel_for))/(max(chunk_behaviourData.vel_for) - min(chunk_behaviourData.vel_for));;
%     
    [B_fRate_temp,FitInfo_fRate] = lasso(X_fRate, y_fRate,'NumLambda',2, 'Alpha', 0.5) ;
    [B_Vm_temp,FitInfo_Vm] = lasso(X_Vm, y_Vm, 'NumLambda',2, 'Alpha',0.5) ;
        
    B_fRate_all(:,count) = B_fRate_temp(:,1); 
    B_Vm_all(:,count) = B_Vm_temp(:,1);

    yCalc_fRate = FitInfo_fRate.Intercept(1) + X_fRate*B_fRate_all(:,count);
    
    yCalc_Vm = FitInfo_Vm.Intercept(1) + X_Vm*B_Vm_all(:,count);
    

    figure();
    subplot(2,1,1)
    normplot(y_Vm - yCalc_Vm)
    subplot(2,1,2)
    normplot(y_fRate - yCalc_fRate)

    R2_Vm_chunk =  1 - sum((y_Vm - yCalc_Vm).^2)/sum((y_Vm - mean(y_Vm)).^2);
    R2_fRate_chunk = 1 - sum((y_fRate - yCalc_fRate).^2)/sum((y_fRate - mean(y_fRate)).^2);

    R2_Vm(count) = R2_Vm_chunk;
    R2_fRate(count) = R2_fRate_chunk;

    count = count + 1; 
    
%     activity_values = chunk_trialData.smooth_Vm';
%     x_values = chunk_behaviourData.angle';
%     y_values = chunk_behaviourData.vel_for';

    

%     x_edges = [-180:10:180]; %
%     y_edges = [-2:.5:10];
%     
%     [N_Vm, heatmap_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);
%     
%     figure(); clf;
%            
%             imagesc(flip(heatmap_Vm))
%             ylabel('Vf mm/s')
%             xt = get(gca, 'XTick');
%             xtnew = linspace(min(xt), max(xt), 7);                             
%             xtlbl = linspace(-135, 165, numel(xtnew));                  
%             set(gca, 'XTick',xtnew, 'XTickLabel',xtlbl)
%             xlabel('angle')
%             colorbar
%             title('Vm')
%             yt = get(gca, 'YTick');
%             ytnew = linspace(min(yt),max(yt),6);
%             ytlbl = linspace(8.5, -1.5, numel(ytnew));
%             set(gca, 'YTick',ytnew, 'YTickLabel',ytlbl)
%             set(gcf,'color','w')
        
end   
toc 



figure(3);clf;
    subplot(2,1,1)
        plot(R2_Vm)
        text(count-1,max(R2_Vm),num2str(mean(R2_Vm)))
    subplot(2,1,2)
        plot(R2_fRate)
        text(count-1,max(R2_fRate),num2str(mean(R2_fRate)))
figure(4);clf;
%     subplot(2,1,1)
    plot(tt,cosBasis*B_Vm_all(2:end,:),'Color',[0.75 0.75 0.75])
    hold on
    plot(tt,mean(cosBasis*B_Vm_all(2:end,:),2),'k')
    box off
    set(gcf,'color',[1 1 1]); 
    
% subplot(2,1,2)
%     plot(tt,cosBasis*B_fRate_all(2:end,:),'Color',[0.75 0.75 0.75])
%     hold on
%     plot(tt,mean(cosBasis*B_fRate_all(2:end,:),2),'k')
    
figure(5);clf;
    subplot(2,1,1)
    plot(B_Vm_all(1,:))
    subplot(2,1,2)
    plot(B_fRate_all(1,:))

