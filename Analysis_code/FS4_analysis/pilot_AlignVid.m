% Sync ephys data & fictrac video - want to know what ephys data point each
% fictrac frame is associated with 

% 1. Align video & ephys, video starts just before daq acquires
%   frames are discrete obj at lower sample rate - lets downsample ephys to
%   match fictrac rate
%   Then compare video signals (DAT) to DAQ signals to align them
%   Result: each downsampled ephys data point should be matched to a video
%   frame



saveVid = 1; 
% rootDir = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\FS4_exps\20230922\fly_1\cell_1\CL_trial_1';
% video = 'fictrac-raw-20230922_120639.avi';
% dat = 'fictrac-2023 0922_120639.dat';
% frameLog = 'fictrac-vidLogFrames-20230922_120639.txt'; 

%rootDir = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\ephys\newData\20231020\fly_1\cell_1\CL_trial_1';

files = dir(fullfile(rootDir,'fictracData')); 
for f = 1:length(files)
    file = files(f).name;
    if contains(file,'.dat')
        dat = file; 
    elseif contains(file,'raw') && contains(file,'.avi')
        video = file; 
    elseif contains(file,'vidLogFrames')
        frameLog = file;
    end  
end

% video = 'fictrac-raw-20231019_164939.avi';
% dat = 'fictrac-20231019_164939.dat';
% frameLog = 'fictrac-vidLogFrames-20231019_164939.txt'; 


DAT = fullfile(rootDir,'fictracData',dat);
[ftD] = readFictracCSV(DAT);

load(fullfile(rootDir,'trialData.mat'));
load(fullfile(rootDir,'trialMeta.mat'));
load(fullfile(rootDir,'pro_trialData.mat'));
tData = processed_trialData{1};
load(fullfile(rootDir,'behaviorData.mat'));
load(fullfile(rootDir,'rawData.mat'));
ephysSettings

temp_for = rawData( :, settings.raw.x );
% transfrom ficTrac signal into radians  
posRadians = temp_for .* 2 .* pi ./ 10; 
% upwrap position signal
unwrappedPos = unwrap( posRadians );
[disp_for,frameTimes_daq] = resample_new(unwrappedPos,60,(settings.sampRate));

% I think I'll have to do this manually, align signals & finddelay give
% wrong value
figure(); 
plot(disp_for)
hold on
plot(ftD.intX)

DAQvalue = 0; %input('DAQ value (blue): ');
DATvalue = 1; %input('DAT value (red): ');

sigDiff = DATvalue - DAQvalue; 

if sigDiff > 0 
    datFor = ftD.intX(sigDiff:end);
    endDiff = length(datFor) - length(disp_for); 
    datFor = datFor(1:end-endDiff);
else
    datFor = ftD.intX;
    endDiff = length(datFor) - length(disp_for); 
    datFor = datFor(1:end-endDiff);
end
figure(); 
plot(disp_for)
hold on
plot(datFor)

% Alignment isn't perfect, likely due to fictrac real sample rate being
% between 59-60 Hz, rough alignment ok for now (look @ what I did in flyg
% for more precise alignment)

frameLog = readmatrix(fullfile(rootDir,'fictracData',frameLog));
allFrames = ftD.frameCounter; 
frameTimes = (ftD.altTimestamp - ftD.altTimestamp(1))/1e3; 

frameLog(:,2) = zeros(size(frameLog));
temp =1:length(frameLog); 

for f = 1:length(frameLog)
    frame = frameLog(f,1);
    try
        fTime = frameTimes(allFrames == frame);
        frameLog(f,2) = fTime;
        idx = find(allFrames == frame); 
    catch
        fTime = frameTimes(idx) + 1/60;
        frameLog(f,2) = fTime;
    end
    
end

if sigDiff > 0 
    datTimes= frameTimes(sigDiff:end); 
    datTimes = datTimes(1:end-endDiff); 
    datTimeEnd = datTimes(end); 
    datTimeStart = datTimes(1); 
else
    datTimes= frameTimes(1:end); 
    datTimes = datTimes(1:end-endDiff); 
    datTimeEnd = datTimes(end); 
    datTimeStart = datTimes(1); 
end
frames = frameLog(frameLog(:,2) > datTimeStart & frameLog(:,2) < datTimeEnd,1);
times = frameLog(frameLog(:,2) > datTimeStart & frameLog(:,2) < datTimeEnd,2);
vidFrame_indices = temp(frameLog(:,2) > datTimeStart & frameLog(:,2) < datTimeEnd)';

vidFrames = table(frames, times, vidFrame_indices);

mov = VideoReader(fullfile(rootDir,'fictracData',video));

trialMeta.notes

%%

try
    mov_cut = read(mov, [vidFrames.vidFrame_indices(1) vidFrames.vidFrame_indices(end)]);
catch
    %mov_cut = read(mov, [vidFrames.vidFrame_indices(1) size(vidFrames,1)]);
    mov_cut = read(mov);
    mov_cut = mov_cut(:,:,:,vidFrames.vidFrame_indices(1):end); 
end


new_mov = VideoWriter('mov_cut','Motion JPEG AVI');
new_mov.FrameRate = 120; 

figure('Position',[0,200,600,900]);
set(gca,"NextPlot","replacechildren")
ax1 = subplot(6,1,1:4); 
ax2 = subplot(6,1,5:6); 

nFrames = size(mov_cut,4);
nDataPoints = length(tData.scaledOutput);
step = round(nDataPoints/nFrames); 
index = 1:step:nDataPoints; 

% Display the first frame in the top subplot
vidFrame = squeeze(mov_cut(:,:,:,1));
image(vidFrame, 'Parent', ax1);
ax1.Visible = 'off';

i=2;
h = plot(ax2,tData.time(1:index(i)), tData.scaledOutput(1:index(i)),'k');
%hold on;
%dot = plot(ax2,tData.time(index(i)),tData.scaledOutput(index(i)),'ro');
%ax2.XLim = [tData.time(1) tData.time(end)];
ax2.YLim = [min(tData.scaledOutput) max(tData.scaledOutput)];
line = xline(ax2, tData.time(index(2)),'r');

for v = 1:nFrames
    pause(1/new_mov.FrameRate);
    
    vidFrame = squeeze(mov_cut(:,:,:,v)); 
    image(vidFrame,'Parent',ax1);
    ax1.Visible = 'off';
    
    i = i+1;
    window = 10*1000; 
    

    if index(i) - window/2 < 1
        si = 1;
    else
        si = index(i) - window/2;
    end


    if index(i) + window/2 > length(tData.time)
        ei = length(tData.time);
    else
        ei = index(i) + window/2;
    end
   
    
    set(h,'YData',tData.scaledOutput(si:ei),'XData',tData.time(si:ei));   
    delete(line);
    line = xline(ax2, tData.time(index(i)),'r');
    set(line);

    %set(dot, 'YData',tData.scaledOutput(index(i)),'XData',tData.time(index(i)))
    ax2.XLim = [tData.time(si) tData.time(ei)]; 
    ax2.YLim = [-60 -10]; 
end
%%
tic

timeStart = 539;%min(vidFrames.times);
startIdx = find(vidFrames.times >= timeStart,1,'first');

timeEnd = 590;%max(vidFrames.times); 
endIdx = find(vidFrames.times <= timeEnd,1,'last');

if saveVid
    mov_cut = read(mov, [vidFrames.vidFrame_indices(startIdx) vidFrames.vidFrame_indices(endIdx)]);
    new_mov = VideoWriter('mov_cut','Motion JPEG AVI');
    new_mov.FrameRate = 120; 


    vid = VideoWriter(fullfile(rootDir,'figures','wholeTrial_movie_Q.avi'));
    open(vid);

    f = figure('Position',[0,200,600,900]);%,'visible','on');
    ax1 = subplot(6,1,1:4); 
    ax2 = subplot(6,1,5:6); 

    nFrames = size(mov_cut,4);
    
    startIdx = find(tData.time >= timeStart,1,'first');
    endIdx = find(tData.time <= timeEnd,1,'last');
    
    activity = tData.scaledOutput(startIdx:endIdx);
    time = tData.time(startIdx:endIdx);
    
    nDataPoints = length(activity);
    step = round(nDataPoints/nFrames); 
    index = 1:step:nDataPoints; 

    % Display the first frame in the top subplot
    vidFrame = squeeze(mov_cut(:,:,:,1));
    image(vidFrame, 'Parent', ax1);
    ax1.Visible = 'off';

    i=2;
    h = plot(ax2,time(1:index(i)), activity(1:index(i)),'k');
    %hold on;
    %dot = plot(ax2,tData.time(index(i)),tData.scaledOutput(index(i)),'ro');
    %ax2.XLim = [tData.time(1) tData.time(end)];
    %ax2.YLim = [min(activity) max(activity)];
    ax2.YLim = [-60 -20];
    line = xline(ax2, time(index(2)),'r');
    
    frameArray = repmat(getframe(f).cdata,[1,1,1,nFrames]);
    for v = 1:nFrames
        %pause(1/new_mov.FrameRate);

        vidFrame = squeeze(mov_cut(:,:,:,v)); 
        image(vidFrame,'Parent',ax1);
        ax1.Visible = 'off';

        i = i+1;
        window = 10*1000; 


        try
            if index(i) - window/2 < 1
                si = 1;
            else
                si = index(i) - window/2;
            end


            if index(i) + window/2 > length(time)
                ei = length(time);
            else
                ei = index(i) + window/2;
            end


            set(h,'YData',activity(si:ei),'XData',time(si:ei));   
            delete(line);
            line = xline(ax2,time(index(i)),'r');
            set(line);

            %set(dot, 'YData',tData.scaledOutput(index(i)),'XData',tData.time(index(i)))
            ax2.XLim = [time(si) time(ei)];
            frameArray(:,:,:,v) = getframe(f).cdata;
        catch
            disp('End of Video')
        end
        
    end
    writeVideo(vid,frameArray); 
    close(vid);
    toc
end





