%% for a given trial this plots the normalized Vm & fRate of the cell to give a rough est of pref heading

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

edges = [-180:10:180];
    [centers_Vm, mean_bin_Vm] = create_binned_mean(angle, Vm, edges);
    [centers_fRate, mean_bin_fRate] = create_binned_mean(angle, fRate_sec, edges);
    figure(1); clf;
    plot(centers_Vm, normalize(mean_bin_Vm),'k')
    hold on 
    plot(centers_fRate, normalize(mean_bin_fRate),'r')
    ylabel('activity')
    xlabel('bar angle')
% look at plot and decide what to call the pref heading & the range of
% headings to take around it
%%
prefHead = 50;
range = 60;

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




% for all time points that fall within this heading range extract out Vf & Vs & activity values
VsinHead = []; 
VfinHead = [];
VminHead = []; 
fRateinHead = []; 
[index, ~, ~] = find(angle > prefHead-range/2 & angle < prefHead+range/2);
[~,~,headingValues] = find(angle(angle > prefHead-range/2 & angle < prefHead+range/2));


%Extract timepoints within heading range where Vs in either only - or +
VsinHead = Vs(index); 
VfinHead = Vf(index);
VminHead = Vm(index); 
fRateinHead = fRate_sec(index);
% Project vel w/ matching side directionality to projection vector
% + = right, - = left
rVs = VsinHead(VsinHead > 0);
rVf = VfinHead(VsinHead > 0);
rVm = VminHead(VsinHead > 0);
rfRate = fRateinHead(VsinHead > 0); 

lVs = VsinHead(VsinHead < 0);
lVf = VfinHead(VsinHead < 0);
lVm = VminHead(VsinHead < 0);
lfRate = fRateinHead(VsinHead < 0); 

%velDot(:,1) = right, velDot(:,2) = left
velDataR = [];
velDataL = [];
for d = 1:2
    if d == 1
       velDataR(:,1) = rVs;
       %velDataR(:,1) = 0;
       velDataR(:,2) = rVf;
       %velDataR(:,2) = 0;
       deg = 31; 
       velDotR = zeros(length(rVs),1);
        for i = 1:length(rVs) 
            velDotR(i,1) = dot(velDataR(i,[1 2]), [sind(deg) cosd(deg)]); 
        end
        
        edges = [min(velDotR(:,1)):0.1:max(velDotR(:,1))]; 
        [N, edges, bin] = histcounts(velDotR(:,1), edges);
        temp = accumarray(bin+1, rfRate, [length(edges) 1]);
        mean_binR = bsxfun(@rdivide, temp(2:end), N');
        centersR = edges(1:end-1)+diff(edges)/2;
        
        
    else
       velDataL(:,1) = lVs; 
       %velDataL(:,1) = 0;
       velDataL(:,2) = lVf;
       %velDataL(:,2) = 0;
       deg = -31;
       velDotL = zeros(length(lVs),1);
       
        for i = 1:length(lVs) 
            velDotL(i,1) = dot(velDataL(i,[1 2]), [sind(deg) cosd(deg)]); 
        end
        
        edges = [min(velDotL(:,1)):0.1:max(velDotL(:,1))]; 
        [N, edges, bin] = histcounts(velDotL(:,1), edges);
        temp = accumarray(bin+1, lfRate, [length(edges) 1]);
        mean_binL = bsxfun(@rdivide, temp(2:end), N');
        centersL = (edges(1:end-1)+diff(edges)/2);

    end
end

figure(7); clf;
plot(centersL,mean_binL,'k')
hold on
plot(centersR,mean_binR,'r')
%% Project vel onto +/- 31 regardless of sign  

prefHead = 0;
range = 360;

edges = [prefHead-range/2:5:prefHead+range/2];

[centers_Vm, mean_bin_Vm] = create_binned_mean(angle, Vm, edges);
[centers_fRate, mean_bin_fRate] = create_binned_mean(angle, fRate_sec, edges);
figure(1); clf;
plot(centers_Vm, normalize(mean_bin_Vm),'k')
hold on 
plot(centers_fRate, normalize(mean_bin_fRate),'r')
ylabel('activity')
xlabel('bar angle')

VsinHead = []; 
VfinHead = [];
VminHead = []; 
fRateinHead = []; 
[index, ~, ~] = find(angle > prefHead-range/2 & angle < prefHead+range/2);
[~,~,headingValues] = find(angle(angle > prefHead-range/2 & angle < prefHead+range/2));


%Extrat timepoints within heading range where Vs in either only - or +
VsinHead = Vs(index); 
VfinHead = Vf(index);
VminHead = Vm(index); 
fRateinHead = fRate_sec(index);


velDataR = [];
velDataL = [];
bin = [];
for d = 1:2
    if d == 1
       velDataR(:,1) = VsinHead;
       %velDataR(:,1) = 0;
       velDataR(:,2) = VfinHead;
       %velDataR(:,2) = 0;
       deg = 31; 
       velDotR = zeros(length(VsinHead),1);
        for i = 1:length(VsinHead) 
            velDotR(i,1) = dot(velDataR(i,[1 2]), [sind(deg) cosd(deg)]); 
        end
        
        edges = [min(velDotR(:,1)):0.1:max(velDotR(:,1))]; 
        [N, edges, bin] = histcounts(velDotR(:,1), edges);
        temp = accumarray(bin+1, VminHead, [length(edges) 1]);
        mean_binR = bsxfun(@rdivide, temp(2:end), N');
        centersR = edges(1:end-1)+diff(edges)/2;
        
    else
       velDataL(:,1) = VsinHead; 
       %velDataL(:,1) = 0;
       velDataL(:,2) = VfinHead;
       %velDataL(:,2) = 0;
       deg = -31;
       velDotL = zeros(length(VsinHead),1);
       
        for i = 1:length(VsinHead) 
            velDotL(i,1) = dot(velDataL(i,[1 2]), [sind(deg) cosd(deg)]); 
        end
        
        edges = [min(velDotL(:,1)):0.1:max(velDotL(:,1))]; 
        [N, edges, bin] = histcounts(velDotL(:,1), edges);
        temp = accumarray(bin+1, VminHead, [length(edges) 1]);
        mean_binL = bsxfun(@rdivide, temp(2:end), N');
        centersL = (edges(1:end-1)+diff(edges)/2);
        
    end
end

figure(7); clf;
plot(centersL,mean_binL,'k')
hold on
plot(centersR,mean_binR,'r')

velDataR = [];
velDataL = [];
bin = [];
for d = 1:2
    if d == 1
       velDataR(:,1) = VsinHead;
       %velDataR(:,1) = 0;
       velDataR(:,2) = VfinHead;
       %velDataR(:,2) = 0;
       deg = 31; 
       velDotR = zeros(length(VsinHead),1);
        for i = 1:length(VsinHead) 
            velDotR(i,1) = dot(velDataR(i,[1 2]), [sind(deg) cosd(deg)]); 
        end
        
        edges = [min(velDotR(:,1)):0.1:max(velDotR(:,1))]; 
        [N, edges, bin] = histcounts(velDotR(:,1), edges);
        temp = accumarray(bin+1, fRateinHead, [length(edges) 1]);
        mean_binR = bsxfun(@rdivide, temp(2:end), N');
        centersR = edges(1:end-1)+diff(edges)/2;
        
    else
       velDataL(:,1) = VsinHead; 
       %velDataL(:,1) = 0;
       velDataL(:,2) = VfinHead;
       %velDataL(:,2) = 0;
       deg = -31;
       velDotL = zeros(length(VsinHead),1);
       
        for i = 1:length(VsinHead) 
            velDotL(i,1) = dot(velDataL(i,[1 2]), [sind(deg) cosd(deg)]); 
        end
        
        edges = [min(velDotL(:,1)):0.1:max(velDotL(:,1))]; 
        [N, edges, bin] = histcounts(velDotL(:,1), edges);
        temp = accumarray(bin+1, fRateinHead, [length(edges) 1]);
        mean_binL = bsxfun(@rdivide, temp(2:end), N');
        centersL = (edges(1:end-1)+diff(edges)/2);
        
    end
end


figure(8); clf;
plot(centersL,mean_binL,'k')
hold on
plot(centersR,mean_binR,'r')

%% Bin velocity vectors based on the standard deviation of Vs
clear

rootPath = '/Users/ewesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys'; %'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'; %change dep on comp
date = input('Date? ','s');
cell_num = input('Cell? ','s');
cell_num = strcat('cell_',cell_num);
trial = input('Trial? ','s');
trial = strcat('trial_',trial);
fileName = fullfile(rootPath,date,cell_num,trial);
figureName = strcat(date,' ',cell_num,' ',trial);
cd(fileName)

%load('rawData.mat');
load('pro_trialData.mat');
load('pro_behaviourData.mat');

cd('/Users/ewesteinde/Documents/EphysCode');
%%
%deg = ; 
count = 1;
for deg = -45:1:45

prefHead = 0; %[45, 135, -135, -45];
range = 360;
step = 0.5;
%for prefHead = headings
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


    % edges =[prefHead-range/2:5:prefHead+range/2];
    % 
    % [centers_Vm, mean_bin_Vm] = create_binned_mean(angle, Vm, edges);
    % [centers_fRate, mean_bin_fRate] = create_binned_mean(angle, fRate_sec, edges);
    % figure(1); clf;
    % plot(centers_Vm, normalize(mean_bin_Vm),'k')
    % hold on 
    % plot(centers_fRate, normalize(mean_bin_fRate),'r')
    % ylabel('activity')
    % xlabel('bar angle')

    % % sanity check, does activity at all seem to increase with increasing vel?
    % x_edges = [prefHead-range/2:1:prefHead+range/2];
    % y_edges = [-2:1:8];
    % 
    % [N_fR, heatmap_array_FR, x_centers_FR, y_centers_FR] = create_activity_heatmap(angle, Vf, fRate_sec, x_edges, y_edges);
    % [N_Vm, heatmap_array_Vm, x_centers_Vm, y_centers_Vm] = create_activity_heatmap(angle, Vf, Vm, x_edges, y_edges);
    % 
    % figure(1);clf;
    %     imagesc(flip(heatmap_array_FR))
    %     colorbar
    % figure(2);clf;
    %     imagesc(flip(heatmap_array_Vm))
    %     colorbar   


    VsinHead = []; 
    VfinHead = [];
    VminHead = []; % angle > prefHead-range/2 & angle < prefHead+range/2
    fRateinHead = []; %angle <= -150 | angle >= 140
    [index, ~, ~] = find(angle > prefHead-range/2 & angle < prefHead+range/2);
    [~,~,headingValues] = find(angle(angle > prefHead-range/2 & angle < prefHead+range/2));


    %Extrat timepoints within heading range where Vs in either only - or +
    VsinHead = Vs(index); 
    VfinHead = Vf(index);
    VminHead = Vm(index); 
    fRateinHead = fRate_sec(index);
    % Project vel w/ matching side directionality to projection vector
    % + = right, - = left
    S = std(VsinHead); 

    hVs = VsinHead(VsinHead > S);
    hVf = VfinHead(VsinHead > S);
    hVm = VminHead(VsinHead > S);
    hfRate = fRateinHead(VsinHead > S);

    oVs = VsinHead(VsinHead < -S);
    oVf = VfinHead(VsinHead < -S);
    oVm = VminHead(VsinHead < -S);
    ofRate = fRateinHead(VsinHead < -S);

    zVs = VsinHead(VsinHead > -S & VsinHead < S);
    zVf = VfinHead(VsinHead > -S & VsinHead < S);
    zVm = VminHead(VsinHead > -S & VsinHead < S);
    zfRate = fRateinHead(VsinHead > -S & VsinHead < S);

    velDatah = [];
    velDatao = [];
    velDataz = []; 
    for d = 1:3
        if d == 1
           velDatah(:,1) = hVs;
           velDatah(:,2) = hVf;
           velDoth = zeros(length(hVs),1);
            for i = 1:length(hVs) 
                velDoth(i,1) = dot(velDatah(i,[1 2]), [sind(deg) cosd(deg)]); 
            end

            edges = [min(velDoth(:,1)):step:max(velDoth(:,1))]; 
            [N, edges, bin] = histcounts(velDoth(:,1), edges);
            temp = accumarray(bin+1, hfRate, [length(edges) 1]);
            mean_binh = bsxfun(@rdivide, temp(2:end), N');
            centersh = edges(1:end-1)+diff(edges)/2;


        elseif d == 2
           velDataz(:,1) = zVs; 
           velDataz(:,2) = zVf;
           velDotz = zeros(length(zVs),1);

            for i = 1:length(zVs) 
                velDotz(i,1) = dot(velDataz(i,[1 2]), [sind(deg) cosd(deg)]); 
            end
            bin = [];
            edges = [min(velDotz(:,1)):step:max(velDotz(:,1))]; 
            [N, edges, bin] = histcounts(velDotz(:,1), edges);
            temp = accumarray(bin+1, zfRate, [length(edges) 1]);
            mean_binz = bsxfun(@rdivide, temp(2:end), N');
            centersz = (edges(1:end-1)+diff(edges)/2);

        else
           velDatao(:,1) = oVs; 
           velDatao(:,2) = oVf;
           velDoto = zeros(length(oVs),1);

            for i = 1:length(oVs) 
                velDoto(i,1) = dot(velDatao(i,[1 2]), [sind(deg) cosd(deg)]); 
            end
            bin = [];
            edges = [min(velDoto(:,1)):step:max(velDoto(:,1))]; 
            [N, edges, bin] = histcounts(velDoto(:,1), edges);
            temp = accumarray(bin+1, ofRate, [length(edges) 1]);
            mean_bino = bsxfun(@rdivide, temp(2:end), N');
            centerso = (edges(1:end-1)+diff(edges)/2);

        end
    end
%end

% figure(11); clf;
% plot(centersh,mean_binh,'r')
% title('111420 Cell1 Trial 1 allheadings')
% hold on
% plot(centersz,mean_binz,'k')
% plot(centerso,mean_bino,'b')

h = corrcoef(hVm, velDoth);
z = corrcoef(zVm, velDotz);
o = corrcoef(oVm, velDoto);
sumVectorAngle(count,:) = [deg h(2,1) z(2,1) o(2,1)];
count = count+1; 
end

%%
count = 1;
for deg = -90:1:0
    prefHead = 90; %[45, 135, -135, -45];
    range = 60;
    step = 0.5;
    %for prefHead = headings
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
    VminHead = []; % angle > prefHead-range/2 & angle < prefHead+range/2
    fRateinHead = []; %angle <= -150 | angle >= 140
    [index, ~, ~] = find(angle > prefHead-range/2 & angle < wrapTo180(prefHead+range/2));
    [~,~,headingValues] = find(angle(angle > prefHead-range/2 & angle < prefHead+range/2));



    VsinHead = Vs(index); 
    VfinHead = Vf(index);
    VminHead = Vm(index); 
    fRateinHead = fRate_sec(index);

    rVs = VsinHead;
    rVf = VfinHead;
    rVm = VminHead;
    rfRate = fRateinHead;

    velData = [];
    velData(:,1) = rVs;

    velData(:,2) = rVf;

    for a = 1:2
        if a == 1
           velDot = zeros(length(rVs),1);
            for i = 1:length(rVs) 
                velDot(i,1) = dot(velData(i,[1 2]), [sind(deg) cosd(deg)]); 
            end

            f = corrcoef(rfRate, velDot);
            sumVectorAngle2(count,[1 2]) = [deg f(2,1)];


            edges = [min(velDot(:,1)):0.1:max(velDot(:,1))]; 
            [N, edges, bin] = histcounts(velDot(:,1), edges);
            temp = accumarray(bin+1, rfRate, [length(edges) 1]);
            mean_bin = bsxfun(@rdivide, temp(2:end), N');
            centers = edges(1:end-1)+diff(edges)/2;
%             figure(1);
%             hold on
%             plot(centers(centers >= 0), mean_bin(centers >= 0)')
            c = polyfit(centers(centers >= 0), mean_bin(centers >= 0)',1); 
            y_est = polyval(c,centers(centers >= 0));
%             plot(centers(centers >= 0),y_est,'r--')
            sumVectorAngle2(count,[3 4]) = [c(1) c(2)]; 
        else
            velDot = zeros(length(rVs),1);
            for i = 1:length(rVs) 
                velDot(i,1) = dot(velData(i,[1 2]), [sind(deg) cosd(deg)]); 
            end
            
            v = corrcoef(rVm, velDot);
            sumVectorAngle2(count,5) = v(2,1);


            edges = [min(velDot(:,1)):0.1:max(velDot(:,1))]; 
            [N, edges, bin] = histcounts(velDot(:,1), edges);
            temp = accumarray(bin+1, rVm, [length(edges) 1]);
            mean_bin = bsxfun(@rdivide, temp(2:end), N');
            centers = edges(1:end-1)+diff(edges)/2;
%             figure(2);
%             hold on
%             plot(centers(centers >= 0), mean_bin(centers >= 0)')
            c = polyfit(centers(centers >= 0), mean_bin(centers >= 0)',1); 
            y_est = polyval(c,centers(centers >= 0));
%             plot(centers(centers >= 0),y_est,'r--')
            sumVectorAngle2(count,[6 7]) = [c(1) c(2)]; 
        end
    end
   
        
        
count = count+1;
end
        
disp(find(sumVectorAngle2(:,3) == max(sumVectorAngle2(:,3))))
disp(find(sumVectorAngle2(:,6) == max(sumVectorAngle2(:,6))))


