function [processed_trialData] =  calc_fRate(trialData_all, downsample_Hz, folder)

ephysSettings

% default spike parameters to start with
minPeakHeight = 5;
posThres = 0.5;
negThres = 0.5;

if iscell(trialData_all)
    trials = length(trialData_all);
else
    trials = 1;
end

for t = 1:trials
%% Spike Detection
    % spike detection and Vm filtering optomized for PFL2/3 oscillatory
    % activity 5/20/21, might need to be adjusted for other cell types
    if iscell(trialData_all)
        trialData = trialData_all{t}; 
    else
        trialData = trialData_all; 
    end

     scaledOutput = trialData.scaledOutput;
    [scaledOutput,time] = resample_new(scaledOutput, downsample_Hz, settings.sampRate);
    
    %correct for leak junction potential - added 3/26/21 any cells before
    %this will still need this done
    scaledOutput = scaledOutput-13;
    medianfilteredOutput = medfilt1(scaledOutput, 500); % 500ms time window
    medianfilteredOutput_Vm = medfilt1(scaledOutput, 50);
    subtractedOutput = scaledOutput - medianfilteredOutput;
    smoothVm = smoothdata(medianfilteredOutput_Vm,'loess',20);
     
%     spikes = zeros(size(subtractedOutput));
%     [pks, locs] = findpeaks(subtractedOutput,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDist);
%     spikes(locs) = 1;

    
    
    %% new spike detection 5/20/21
    satisfied = 0; 
    while satisfied == 0 
        %take derivative of the subtracted output
         Vm_derivative = gradient( subtractedOutput ) ; 
        % find all points where the derivative passes a threshold of 0.4
        [~, locs_pos] = findpeaks(Vm_derivative,'MinPeakHeight', posThres);
        % find all points where the derivative passes a threshold of -0.8
        [~, locs_neg] = findpeaks(Vm_derivative .* -1 ,'MinPeakHeight', negThres,'MinPeakDistance', 10);
        % find all points where Vm passes a threshold of 2
        [~, locs_amp] = findpeaks(subtractedOutput,'MinPeakHeight', minPeakHeight);
            % threshold values might need to be adjusted for every
            % experiment
        % If the time difference between a positive deflection and a
        % following negative deflection is 5ms (indexPos-indexNeg <= 5) or less then this is a spike
        % keep positive deflection points where this occurs 
        locs_temp = []; 
        for i = 1:length(locs_neg)
            timePoint = locs_neg(i); 
            tdiff = locs_pos - timePoint;
            spike_start = locs_pos(tdiff >= -10 & tdiff < 0); 
            if length(spike_start) > 1
                spike_start = spike_start(1);
            end
            if ~isempty(spike_start)
                locs_temp = cat(1,locs_temp, spike_start);
            end
        end
        
        locs = [];
        for i = 1:length(locs_temp)
            timePoint = locs_temp(i); 
            tdiff = locs_amp - timePoint; 
            if ~isempty(tdiff(abs(tdiff) <= 10))
                locs = cat(1, locs, timePoint);
            end
        end
        
        spikes = zeros(size(subtractedOutput));
        spikes(locs) = 1;
    %% Spike detection troubleshooting
    figure(8); clf;

                g(1) = subplot(4,1,1);
                plot(time, subtractedOutput, 'k')
                ylabel('sub Output')
                
                g(2) = subplot(4,1,2);
                plot(time, Vm_derivative, 'k')
                ylabel('Vm derivative')
   
                g(3) = subplot(4,1,3);
                plot(time, spikes, 'k')
                ylabel('Spikes')

                g(4) = subplot(4,1,4);
                plot(time, smoothVm, 'k')
                ylabel('Vm')

                linkaxes(g,'x');

    adjust = input('adjust spike detection?','s');

    if strcmp(adjust,'n')
        satisfied = 1;
    else
        minPeakHeight = input('Min peak height? ');
        posThres = input('Positive Threshold? ');
        negThres = input('Negative Threshold? ');
    end
    
    end
    %% Asa method Guassian Kernal to calc continuous firing rate
    
    kernWidth = 20; %ms
    smoothSamples = kernWidth*1; % kernWidth * SampRate (in same units as kernWidth)
    raster = spikes;
    wingSize = 5.5; % Std
    gSamp = [-(wingSize)*smoothSamples:wingSize*smoothSamples];
    %eSamp = [-smoothSamples:smoothSamples];
    %kernal1 = 1/(smoothSamples *
    %sqrt(2*pi)).*exp(-gSamp.^2/(2*smoothSamples.^2)); gaussian kernel 
    kernel = 1/(smoothSamples * sqrt(2)).*exp(-sqrt(2)*abs(gSamp/smoothSamples)); %exponential kernel
    psth = conv(raster,kernel,'same');
    
    disp(sum(kernel)) 

    % multiply by sample rate (in sec) to get correct units
    fRate_sec = psth .* 1000; 
    figure(9); clf;

            h(1) = subplot(4,1,1);
            plot(trialData.time, trialData.scaledOutput,'k')
            ylabel('Vm')

            h(2) = subplot(4,1,2);
            plot(time,psth, 'k')
            ylabel('PSTH')
            
            h(3) = subplot(4,1,3);
            plot(time,fRate_sec, 'k')
            ylabel('fRate_sec')

            h(4) = subplot(4,1,4);
            plot(time,smoothVm, 'k')
            ylabel('Vm, no spikes')

            xlabel('Time (s)')

            linkaxes(h,'x');

%     adjust2 = input('adjust Guassian kernal?','s');
%             
%     if strcmp(adjust2,'y')
%         error('Adjust guassian kernal parameters')
%     end
%     
   
        processed_trialData{t} = table(time, scaledOutput,spikes, subtractedOutput, psth, fRate_sec, smoothVm);

end
    save(fullfile(folder,'pro_trialData.mat'),'processed_trialData')
end