function [pulse_table] = detect_iontoPulses(processed_behaviourData, pulse_window, sampRate)

%% downsampled jump idx
    stim = processed_behaviourData.stim; 

    count = 1; 
    prev_jump = 0; 
    transition = zeros(size(stim)); 

    for idx = 1:length(stim)
        value = stim(idx);
        if count ~= 1 
            if value ~= prev_val && (count - prev_jump > 5) 
                transition(count) = 1;
                prev_jump = count; 
            else
                transition(count) = 0; 
            end

            if ~(isnan(value))
                prev_val = value;
            else
                prev_val = prev_val;
            end
            
            count = count + 1; 
        else
            transition(count) = 0; 
            prev_val = value;
            count = count + 1;
        end
    end
 
    window_datapoints = round(pulse_window*sampRate); 
    pulse_idx = find(transition == 1);
    count = 1; 
    for idx = 1:length(pulse_idx)
        if rem(idx,2) ~= 0 
            pulse_StartStopidx(count,1) = pulse_idx(idx);
        else
            pulse_StartStopidx(count,2) = pulse_idx(idx); 
            count = count + 1; 
        end
    end
    
    pulse_StartStopidx(:,3) = pulse_StartStopidx(:,2) - pulse_StartStopidx(:,1); 

    count = 1; 
    pulse_array = zeros(size(pulse_StartStopidx,1),5); 
    for idx = 1:size(pulse_StartStopidx,1)
        pre_jump = pulse_StartStopidx(idx,1) - window_datapoints;
        post_jump = pulse_StartStopidx(idx,1) + window_datapoints;
        if post_jump > length(processed_behaviourData.stim)
            post_jump = length(processed_behaviourData.stim);
        end
        pulse_array(count,1) = pre_jump; 
        pulse_array(count,2) = pulse_StartStopidx(idx,1); 
        pulse_array(count,3) = pulse_StartStopidx(idx,2); 
        pulse_array(count,4) = post_jump; 
        pulse_array(count,5) = pulse_StartStopidx(idx,3); 
        count = count + 1;
    end
    
    headers = {'windowStart', 'pulseStart','pulseEnd','windowEnd','pulseLength'};
    pulse_table = array2table(pulse_array,'VariableNames',headers); 
    
end