%% Calculate correlation b/w Vm/FR & Vs/Vf
% Use Pearson's correlation --> X & Y should come from a normal
% distribution 
%% load in the data, pick a preferred direction and range & only look at time points that fall within this preferred heading range
 % +/-30 -150 best corr b/w FR/Vm & Vf
prefHead = 0;
range = 360;
step = 0.1;

Vm = [];
fRate_sec = [];
angle = [];
Vf = [];
Vs = [];
for t = 1:length(processed_trialData)
        
    fRate_sec = cat(1,fRate_sec, processed_trialData{t}.fRate_sec(0.5*1000:end-(0.5*1000)));
    Vm = cat(1,Vm, processed_trialData{t}.smooth_Vm(0.5*1000:end-(0.5*1000)));
    angle = cat(1, angle, processed_behaviourData{t}.angle(0.5*1000:end-(0.5*1000)));
    Vf = cat(1, Vf, processed_behaviourData{t}.vel_for(0.5*1000:end-(0.5*1000)));
    Vs = cat(1, Vs, processed_behaviourData{t}.vel_side(0.5*1000:end-(0.5*1000)));
end

%correct for leak junction potential 
Vm = Vm -13; 

figure(8);clf;
plot(Vf)
hold on; 
%plot(normalize(Vm))
plot(fRate_sec)
% 
edges =[prefHead-range/2:5:prefHead+range/2];

[centers_Vm, mean_bin_Vm] = create_binned_mean(angle, Vm, edges);
[centers_fRate, mean_bin_fRate] = create_binned_mean(angle, fRate_sec, edges);
figure(1); clf;
plot(centers_Vm, normalize(mean_bin_Vm),'k')
hold on 
plot(centers_fRate, normalize(mean_bin_fRate),'r')
ylabel('activity')
xlabel('bar angle')

% sanity check, does activity at all seem to increase with increasing vel?
x_edges = [prefHead-range/2:1:prefHead+range/2];
y_edges = [-2:1:8];

[N_fR, heatmap_array_FR, x_centers_FR, y_centers_FR] = create_activity_heatmap(angle, Vf, fRate_sec, x_edges, y_edges);
[N_Vm, heatmap_array_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(angle, Vf, Vm, x_edges, y_edges);

figure(1);clf;
    imagesc(flip(heatmap_array_FR))
    colorbar
figure(2);clf;
    imagesc(flip(heatmap_array_Vm))
    colorbar   


VsinHead = []; 
VfinHead = [];
VminHead = []; 
fRateinHead = []; 
[index, ~, ~] = find(angle > prefHead-range/2 & angle < prefHead+range/2);
[~,~,headingValues] = find(angle(angle > prefHead-range/2 & angle < prefHead+range/2));


%Extrat timepoints within heading range 
VsinHead = Vs(index); 
VfinHead = Vf(index);
VminHead = Vm(index); 
fRateinHead = fRate_sec(index);
% look at relationship b/w Vm/FR & Vf
A = VminHead;
B = VfinHead;
% extra timepoints w/ + or - Vs 
[corrCof, p, RL, RU] = corrcoef(A,B); 

[c, lags] = xcorr(A,B,1000);
figure(18); plot(lags./1000,c)
%% compare across time points
count = 1; 
nShift = -200:1:200;
VfShift = zeros(length(VfinHead),length(nShift));
for i = nShift %ms
    VfShift(:,count) = circshift(VfinHead,i); 
    if i > 0 
        VfShift(1:i,count) = NaN; 
    elseif i < 0
       VfShift([end + i:end],count) = NaN; 
    end
    count = count + 1; 
end
sCorrP = zeros(length(VfShift(1,:)),4); 
for s = 1:length(VfShift(1,:))
    [corrCof, p, RL, RU] = corrcoef(fRateinHead,VfShift(:,s),'Rows','pairwise'); 
    sCorrP(s,1) = corrCof(2,1);
    sCorrP(s,2) = p(2,1);
    sCorrP(s,3) = RL(2,1);
    sCorrP(s,4) = RU(2,1);
end


%% compare velproj within a set heading range, no time shifts
 % +/-30 -150 best corr b/w FR/Vm & Vf
prefHead = 0;
range = 360;
step = 0.1;

Vm = [];
fRate_sec = [];
angle = [];
Vf = [];
Vs = [];
for t = 1:length(processed_trialData)
        
    fRate_sec = cat(1,fRate_sec, processed_trialData{t}.fRate_sec(0.5*1000:end-(0.5*1000)));
    Vm = cat(1,Vm, processed_trialData{t}.smooth_Vm(0.5*1000:end-(0.5*1000)));
    angle = cat(1, angle, processed_behaviourData{t}.angle(0.5*1000:end-(0.5*1000)));
    Vf = cat(1, Vf, processed_behaviourData{t}.vel_for(0.5*1000:end-(0.5*1000)));
    Vs = cat(1, Vs, processed_behaviourData{t}.vel_side(0.5*1000:end-(0.5*1000)));
end

%correct for leak junction potential 
Vm = Vm -13; 

figure(8);clf;
plot(Vf)
hold on; 
%plot(normalize(Vm))
plot(fRate_sec)

edges =[prefHead-range/2:5:prefHead+range/2];

[centers_Vm, mean_bin_Vm] = create_binned_mean(angle, Vm, edges);
[centers_fRate, mean_bin_fRate] = create_binned_mean(angle, fRate_sec, edges);
figure(1); clf;
plot(centers_Vm, normalize(mean_bin_Vm),'k')
hold on 
plot(centers_fRate, normalize(mean_bin_fRate),'r')
ylabel('activity')
xlabel('bar angle')

% sanity check, does activity at all seem to increase with increasing vel?
x_edges = [prefHead-range/2:1:prefHead+range/2];
y_edges = [-2:1:8];

[N_fR, heatmap_array_FR, x_centers_FR, y_centers_FR] = create_activity_heatmap(angle, Vf, fRate_sec, x_edges, y_edges);
[N_Vm, heatmap_array_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(angle, Vf, Vm, x_edges, y_edges);

figure(1);clf;
    imagesc(flip(heatmap_array_FR))
    colorbar
figure(2);clf;
    imagesc(flip(heatmap_array_Vm))
    colorbar   


VsinHead = []; 
VfinHead = [];
VminHead = []; 
fRateinHead = []; 
[index, ~, ~] = find(angle >= prefHead-range/2 & angle <= prefHead+range/2);
[~,~,headingValues] = find(angle(angle > prefHead-range/2 & angle < prefHead+range/2));


%Extrat timepoints within heading range 
VsinHead = Vs(index); 
VfinHead = Vf(index);
VminHead = Vm(index); 
fRateinHead = fRate_sec(index);

velData = [];
velData(:,1) = VsinHead;
%velDataR(:,1) = 0;
velData(:,2) = VfinHead;
%velDataR(:,2) = 0;
deg = -31; 
velDot = zeros(length(VsinHead),1);
for i = 1:length(VsinHead) 
    velDot(i,1) = dot(velData(i,[1 2]), [sind(deg) cosd(deg)]); 
end
  
[corrCof, p, RL, RU] = corrcoef(fRateinHead,velDot,'Rows','pairwise');
figure(10);clf; plot(velDot); hold on; plot(fRateinHead)

[c, lags] = xcorr(VminHead,velDot,1000,'normalize');
figure(18); plot(lags./1000,c)
%%
find(sCorrP(:,1) == max(sCorrP))


%% compare velproj across times points at all headings
% using the velocity project gives a more consistent correlation shift than
% just using forward velocity, takes into account the influence of Vs
Vm = [];
fRate_sec = [];
angle = [];
Vf = [];
Vs = [];
for t = 1:length(processed_trialData)
        
    fRate_sec = cat(1,fRate_sec, processed_trialData{t}.fRate_sec(0.5*1000:end-(0.5*1000)));
    Vm = cat(1,Vm, processed_trialData{t}.smooth_Vm(0.5*1000:end-(0.5*1000)));
    angle = cat(1, angle, processed_behaviourData{t}.angle(0.5*1000:end-(0.5*1000)));
    Vf = cat(1, Vf, processed_behaviourData{t}.vel_for(0.5*1000:end-(0.5*1000)));
    Vs = cat(1, Vs, processed_behaviourData{t}.vel_side(0.5*1000:end-(0.5*1000)));
end


velData = [];
velData(:,1) = Vs;
%velDataR(:,1) = 0;
velData(:,2) = Vf;
%velDataR(:,2) = 0;
deg = -31; 
velDot = zeros(length(Vf),1);
for i = 1:length(Vs) 
    velDot(i,1) = dot(velData(i,[1 2]), [sind(deg) cosd(deg)]); 
end
        
count = 1; 
nShift = -100:1:100;
VShift = zeros(length(velDot),length(nShift));
for i = nShift %ms
    VShift(:,count) = circshift(velDot,i); 
    if i > 0 
        VShift(1:i,count) = NaN; 
    elseif i < 0
       VShift([end + i:end],count) = NaN; 
    end
    count = count + 1; 
end
sCorrP = zeros(length(VShift(1,:)),4); 
for s = 1:length(VShift(1,:))
    [corrCof, p, RL, RU] = corrcoef(fRate_sec,VShift(:,s),'Rows','pairwise'); 
    sCorrP(s,1) = corrCof(2,1);
    sCorrP(s,2) = p(2,1);
    sCorrP(s,3) = RL(2,1);
    sCorrP(s,4) = RU(2,1);
end
find(sCorrP(:,1) == max(sCorrP))

%% compare velproj across times points within a specific heading range 

prefHead = 30;
range = 120;
step = 0.5;

Vm = [];
fRate_sec = [];
angle = [];
Vf = [];
Vs = [];
for t = 1:length(processed_trialData)
        
    fRate_sec = cat(1,fRate_sec, processed_trialData{t}.fRate_sec(0.5*1000:end-(0.5*1000)));
    Vm = cat(1,Vm, processed_trialData{t}.smooth_Vm(0.5*1000:end-(0.5*1000)));
    angle = cat(1, angle, processed_behaviourData{t}.angle(0.5*1000:end-(0.5*1000)));
    Vf = cat(1, Vf, processed_behaviourData{t}.vel_for(0.5*1000:end-(0.5*1000)));
    Vs = cat(1, Vs, processed_behaviourData{t}.vel_side(0.5*1000:end-(0.5*1000)));
end

%correct for leak junction potential 
Vm = Vm -13; 

VsinHead = []; 
VfinHead = [];
VminHead = []; %angle >= prefHead-range/2 & angle <= prefHead+range/2
fRateinHead = []; % angle <= -150 | angle >= 140
[index, ~, ~] = find(angle >= prefHead-range/2 & angle <= prefHead+range/2); %

%Extract timepoints within heading range 
VsinHead = Vs(index); 
VfinHead = Vf(index);
VminHead = Vm(index); 
fRateinHead = fRate_sec(index);

% identify where indicies are not continuous & sep out trial fragments
count = 1;
start = 0;
frag = {};
shiftPoint = [];
for i = 2:length(index)
    if index(i) - index(i-1) ~= 1
        shiftPoint(1,2) = i-1;
        if count == 1
            frag{count} = index(1:shiftPoint(1,2)); 
            count = count+1;
            shiftPoint(1,1) = i;
        else
            frag{count} = index(shiftPoint(1,1):shiftPoint(1,2));
            count = count+1;
            shiftPoint(1,1) = i;
        end 
    end
    
    if i == length(index)
        if isempty(shiftPoint)
            frag{count} = index(1:i);
        else
            frag{count} = index(shiftPoint(1,1):i);
        end
    end
end
% for each fragment find it's shifted indicies & collect in array
nShift = -100:1:100;
velInd = [];
for g = 1:length(frag)
    count = 1; 
    temp = [];
    for i = nShift %ms
        temp(:,count) = frag{g}(:,1) - i;
        count = count + 1;
    end
    velInd = cat(1,velInd, temp);
end

velData = [];
velData(:,1) = Vs;
%velDataR(:,1) = 0;
velData(:,2) = Vf;
%velDataR(:,2) = 0;
deg = -31; 
velDot = zeros(length(Vf),1);
for i = 1:length(Vs) 
    velDot(i,1) = dot(velData(i,[1 2]), [sind(deg) cosd(deg)]); 
end
velInd(velInd > max(index)) = NaN;
velInd(velInd < min(index)) = NaN;

VShift = zeros(size(velInd)); 
for i = 1:numel(velInd)
    if isnan(velInd(i))
        VShift(i) = NaN;
    else
        VShift(i) = velDot(velInd(i)); 
    end
end

sCorrP = zeros(length(VShift(1,:)),4); 
for s = 1:length(VShift(1,:))
    [corrCof, p, RL, RU] = corrcoef(VminHead,VShift(:,s),'Rows','pairwise'); 
    sCorrP(s,1) = corrCof(2,1);
    sCorrP(s,2) = p(2,1);
    sCorrP(s,3) = RL(2,1);
    sCorrP(s,4) = RU(2,1);
end
iMax = find(sCorrP(:,1) == max(sCorrP))
sCorrP(iMax)
%figure(10);clf; plot(velDot); hold on; plot(fRateinHead)

%%

c = [1 0 0;1 0 0; 1 0 0; 0 1 0; 0 0 1];
sz = 50;
figure(15);clf;
scatter([54,20,37,27,24], [0.528, 0.6253, 0.4756, 0.3855, 0.4062],sz,c,'filled')
xlabel('Delay FR relative to velocity (ms)')
ylabel("Pearson's correlation coefficients")
xlim([0 55])

