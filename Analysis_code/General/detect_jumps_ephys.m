function [jump_array, transition, yChannel] = detect_jumps_ephys(frY, pre_jump_window, post_jump_window, sampRate)

    if size(frY,1) < size(frY,2)
        yChannel = frY'; 
    else 
        yChannel = frY; 
    end 
    
    yChannel_temp = round(yChannel - [1 2 3 4]);
    yChannel_clean = zeros(size(yChannel));
    
    yError = 0; 
    for row = 1:size(yChannel_temp,1)
        yChannel_idx = find(yChannel_temp(row,:) == 0);
        try 
            yChannel_clean(row,1) = yChannel_idx;
        catch
            yError = 1; 
            yChannel_clean(row,1) = 0; 
        end
    end
    
    if yError
        disp('yChannel error')
    end
    
    yChannel = yChannel_clean; 

    count = 1; 
    prev_jump = 0; 
    jcount = 1; 
    transition = zeros(size(yChannel)); 

    for idx = 1:length(yChannel)
        value = yChannel(idx);
        if isnan(value)
            print ('nan approaching')
        end

        if count ~= 1 
            if value ~= prev_val && ~(isnan(value)) && (count - prev_jump > 20) 
                transition(count) = 1;
                prev_jump = count; 
                jumpStep = yChannel_clean(idx+20) - yChannel_clean(idx-20); 
                if jumpStep == 1
                    jumpSize(jcount) = 90; 
                elseif jumpStep == -1
                    jumpSize(jcount) = -90; 
                else
                    jumpSize(jcount) = 180; 
                end
                jcount = jcount + 1; 
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
 
    pre_window_datapoints = round(pre_jump_window*sampRate); 
    post_window_datapoints = round(post_jump_window*sampRate);
    jump_idx = find(transition == 1); 
    jump_array = zeros(length(jump_idx),3); 

    count = 1; 

    for idx = 1:length(jump_idx)
        pre_jump = jump_idx(idx) - pre_window_datapoints;
        post_jump = jump_idx(idx) + post_window_datapoints;
        jump_array(count,1) = pre_jump; 
        jump_array(count,2) = jump_idx(idx); 
        jump_array(count,3) = post_jump; 
        jump_array(count,4) = jumpSize(count);
        count = count + 1;
    end
    
    %Debugging plot 
%     figure()
%     plot(yChannel)
%     hold on
%     plot(transition,'ro')
    
 
    
end