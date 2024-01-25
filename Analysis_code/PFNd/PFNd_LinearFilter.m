%compute a linear filter of PFNd neural activity & translational velocity
%extract timepoints during which the fly's velocity stays above a set
%threshold 

clear

rootPath = '/Users/elenawesteinde/Dropbox (HMS)/Wilson_Lab_Data/ephys'; % 'C:\Users\ewest\Dropbox (HMS)\Wilson_Lab_Data\ephys'; %change dep on comp
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

if isfield(trialMeta, 'notes')
    disp(trialMeta.notes)
end

cd('/Users/elenawesteinde/Documents/EphysCode/Analysis_code'); 
%% extract timepoints when the fly is moving
minVel = 0.5; %mm/s only taking into account side & forward velocities
trial_fragments = cell(1,3); 

for t = 1:length(processed_trialData)
    speed_singleTrial = abs(processed_behaviourData{t}.vel_dot); 
    [index] = find(speed_singleTrial > minVel); 
    count = 1;
    start = 0;
    frag = {};
    shiftPoint = [];
    for i = 2:length(index)
<<<<<<< HEAD
        if index(i) - index(i-1) > 20
=======
        if index(i) - index(i-1) > 200
>>>>>>> 3997c0d47747a7cb3fce5cd2465aa8de4033c30b
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
    
%     for c = 1:length(frag)
%         if frag{c+1}(1) - frag{c}(end) < 500
%             cont_frag{c} = 
    
    trial_fragments{t} = frag; 
end


%% Linear Filter    
% for each continuous segment calculate the linear filter b/w neural
    % activity & the fly's translational velocity
    % r(t) = neural resp over time
    % s(t) = translation vel of ball over time
    % 1. take fast fourier transform of r(t) & s(t) --> S(omega), R(omega) to
    % convert them into the frequency domain 
linearFilterSum = cell(1,3); 
linearFilterAve{t} = cell(1,3); 
for t = 1:length(trial_fragments)
    linearFilterSum{t} = zeros(1,1000);
    count = 0; 
    for i = 1:length(trial_fragments{t})
<<<<<<< HEAD
        
        stimulus = processed_behaviourData{t}.vel_dot(trial_fragments{t}{i}(1):trial_fragments{t}{i}(end));
        response = processed_trialData{t}.fRate_sec(trial_fragments{t}{i}(1):trial_fragments{t}{i}(end));
        
            %uncomment to run over neural activity recording with spikes
    %included
        
%         nActivity = resample(trialData{t}.scaledOutput((trial_fragments{t}{i}(1):trial_fragments{t}{i}(end)), downsample_Hz, settings.sampRate);
%         response = nActivity; 

        stimulusFft = fft(stimulus); 
        responseFft = fft(response);

        figure(1); clf; 
        yyaxis left
        plot(stimulus)
        yyaxis right
        plot(response)
        
        figure(2); clf; 
        yyaxis left
        plot(stimulusFft)
        yyaxis right
        plot(responseFft)
    
    
        % 2. take the autocorrelation of R, self multiplication of R w/ its complex
        % conjugate to retain phase information 
        % self multiplicaiton in freq domain == convolution of R in the time domain
        % 3. multiply R by complexconj(S) to output the relationshio b/w neural resp & the ball
        % 4. divide by R*complexconj(R) to remove similarity that arises from the
        % autocorrelation of R
        
        autocorr = (stimulusFft .* conj(stimulusFft));
        figure(3); clf
        plot(autocorr)
        
        filterFft = (responseFft .* conj(stimulusFft));% ./ (autocorr);
        % take the inverse fast fourier transform of this to return the linear
        % filter
        linearFilter = real(ifft(filterFft));
        
        figure(4); clf
        plot(linearFilter)
        %%
        %how do you sum & average together if all their lengths are
        %dependent on the original chunk length?
        linearFilterSum{i} = linearFilterSum{t} + linearFilter; 
=======
        if length(trial_fragments{t}{1,i}) >= 1000
        %%
            stimulus = processed_behaviourData{t}.vel_dot(trial_fragments{t}{i}(1):trial_fragments{t}{i}(end));
            stimulus = stimulus';
            response = processed_trialData{t}.smooth_Vm(trial_fragments{t}{i}(1):trial_fragments{t}{i}(end));

                %uncomment to run over neural activity recording with spikes
        %included

    %         nActivity = resample(trialData{t}.scaledOutput((trial_fragments{t}{i}(1):trial_fragments{t}{i}(end)), downsample_Hz, settings.sampRate);
    %         response = nActivity; 

            stimulusFft = fft(stimulus); 
            responseFft = fft(response);

            figure(1); clf; 
            yyaxis left
            plot(stimulus)
            yyaxis right
            plot(response)

            figure(2); clf; 
            yyaxis left
            plot(stimulusFft)
            yyaxis right
            plot(responseFft)


            % 2. take the autocorrelation of R, self multiplication of R w/ its complex
            % conjugate to retain phase information 
            % self multiplicaiton in freq domain == convolution of R in the time domain
            % 3. multiply R by complexconj(S) to output the relationshio b/w neural resp & the ball
            % 4. divide by R*complexconj(R) to remove similarity that arises from the
            % autocorrelation of R

            autocorr = (stimulusFft .* conj(stimulusFft));
            figure(3); clf
            plot(autocorr)

            filterFft = (responseFft .* conj(stimulusFft)) %./ (autocorr);
            % take the inverse fast fourier transform of this to return the linear
            % filter
            linearFilter = real(ifft(filterFft));

            figure(4); clf
            plot(linearFilter)
            linearFilterSum{t} = linearFilterSum{t} + linearFilter(1:1000); 
        end
        %%
        %how do you sum & average together if all their lengths are
        %dependent on the original chunk length?
>>>>>>> 3997c0d47747a7cb3fce5cd2465aa8de4033c30b
        count = count + 1; 
    end
    linearFilterAve{t} = linearFilterSum{t}/count; 
    figure(7); clf;
    plot(linearFilterAve{3})
end


% Average together the linear filters of each continous movement segement 