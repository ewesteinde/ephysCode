prefHead = 0; 
step = 0.5;

offset = -31;

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
    
    % correct for leak junction potential 
    Vm = Vm -13;
    
    range = 120;
    head = 120;
    saveBinsVm = cell(1,1);
    saveBinsFR = cell(1,1);
    
    % find indicies of all timepoints that fall within set heading range
    [index, ~, ~] = find(angle > head-range/2 & angle < wrapTo180(head+range/2));
    
    VsinHead = Vs(index);
    VfinHead = Vf(index);
    fRateinHead = fRate_Sec(index);
    VminHead = Vm(index);

% Project Vf & Vs onto Vtpref
Vspref = Vs.*cosd((angle + 90)-wrapTo180(prefHead+offset));
Vfpref = Vf.*cosd(angle-wrapTo180(prefHead+offset));
velDot = Vspref + Vfpref; 

[index, ~, ~] = find(angle > head-range/2 & angle < wrapTo180(head+range/2));
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


% Look at cell of fragments and choose the largest to use as your single
% continuous vector of timepoints
corrArray = [];
count = 0;
for i = 1:length(frag)
    if length(frag{i}) > 1000
        count = count+1;
        index = frag{i}; %set manually

        VsinHead = Vs(index); 
        VfinHead = Vf(index);
        VminHead = Vm(index); 
        fRateinHead = fRate_sec(index);
        velDot = Vspref(index);

            [c, lags] = xcorr(VminHead,velDot,1000,'normalized');
            corrArray(:,count) = c;
    end

end

[c, lags] = xcorr(VminHead,velDot,1000,'normalized');
figure(19); plot(lags./1000,c)
%%
    Vm = [];
    fRate_sec = [];
    angle = [];
    Vf = [];
    Vs = [];
    
    prefHead = 0; 

for t = 1:length(processed_trialData)

    fRate_sec = cat(1,fRate_sec, processed_trialData{t}.fRate_sec(0.5*1000:end-(0.5*1000)));
    Vm = cat(1,Vm, processed_trialData{t}.smooth_Vm(0.5*1000:end-(0.5*1000)));
    angle = cat(1, angle, processed_behaviourData{t}.angle(0.5*1000:end-(0.5*1000)));
    Vf = cat(1, Vf, processed_behaviourData{t}.vel_for(0.5*1000:end-(0.5*1000)));
    Vs = cat(1, Vs, processed_behaviourData{t}.vel_side(0.5*1000:end-(0.5*1000)));
end

%     Vm = Vm(3500:end);
%     fRate_sec = fRate_sec(3500:end);
%     angle = angle(3500:end);
%     Vf = Vf(3500:end);
%     Vs = Vs(3500:end);
%     angle = angle(3500:end);


%correct for leak junction potential 
Vm = Vm-13;

%change based on cell & desired heading range
offset = -31;
head = 0;
range = 360;

% find indicies of all timepoints that fall within set heading range
[index, ~, ~] = find(angle >= head-range/2 | angle <= wrapTo180(head+range/2));

VminHead = Vm(index);
fRateinHead = fRate_sec(index);
VfinHead = Vf(index);
VsinHead = Vs(index);
angleinHead = angle(index);

% clean up later, converts Vf and Vs to Vf and Vs in the direction of the
% preferred heading
Vspref = VsinHead.*cosd(prefHead-angleinHead) + VfinHead.*sind(prefHead-angleinHead);
Vfpref = VfinHead.*cosd(prefHead-angleinHead) + VsinHead.*sind(prefHead-angleinHead);

velData = [];
velData(:,1) = Vspref;
velData(:,2) = Vfpref;

% projects Vf and Vs onto the travel direction vector
velDot = zeros(length(VsinHead),1);
for i = 1:length(VsinHead) 
    velDot(i,1) = dot(velData(i,[1 2]), [sind(offset) cosd(offset)]); 
end

[c, lags] = xcorr(normalize(VminHead),normalize(velDot),1000);
figure(19); plot(lags,c)
lags(c == max(c))
title('120420 Cell1 Trial1 prefAngle 30 range +/-120')
xlabel('time (ms)')

figure(6); clf;
plot(normalize(VminHead),'r')
hold on
plot(velDot,'b')
