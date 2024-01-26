%troubleshooting code calcualtes the linear filter over the entire continous trial rather than during each short segement of continous movement        

<<<<<<< HEAD
        stimulus = processed_behaviourData{t}.vel_dot;
        response = processed_trialData{t}.fRate_sec';
=======
        stimulus = processed_behaviourData{t}.vel_dot';
        response = processed_trialData{t}.smooth_Vm;
>>>>>>> 3997c0d47747a7cb3fce5cd2465aa8de4033c30b
        ephysSettings
        downsample_Hz = 1000;
    %uncomment to run over neural activity recording with spikes
    %included
        
%         nActivity = resample(trialData{t}.scaledOutput, downsample_Hz, settings.sampRate);
%         response = nActivity; 
        
        stimulusFft = fft(stimulus); 
        responseFft = fft(response);

        figure(); clf; 
        yyaxis left
        plot(real(stimulusFft))
        yyaxis right
        plot(real(responseFft))
    
        figure(); clf; 
        yyaxis left
        plot(stimulus)
        yyaxis right
        plot(response)
    
        % 2. take the autocorrelation of R, self multiplication of R w/ its complex
        % conjugate to retain phase information 
        % self multiplicaiton in freq domain == convolution of R in the time domain
        % 3. multiply R by complexconj(S) to output the relationshio b/w neural resp & the ball
        % 4. divide by R*complexconj(R) to remove similarity that arises from the
        % autocorrelation of R
        autocorr = (stimulusFft .* conj(stimulusFft));
        figure(); clf
        plot(autocorr)
        
        filterFft = (responseFft .* conj(stimulusFft));% ./ (autocorr);
        % take the inverse fast fourier transform of this to return the linear
        % filter
        linearFilter = real(ifft(filterFft));
        figure(27); clf
        plot(linearFilter)
        
        %how do you sum & average together if all their lengths are
        %dependent on the original chunk length?
     