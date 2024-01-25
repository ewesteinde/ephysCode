function [no0Vel_frag, no0Vel_idx] = remove0velocity(speedThres, contThres,  processed_behaviourData)

%% extract timepoint indices when the fly is moving
% at some point ask Jenny for her statisitcal version of this to determine
% what is/isn't a start/stop transition 

% speedThres = minimum speed threshold to define when the fly is and isn't moving, only takes into account side & forward velocities
% contThres = threshold for something to count as a continuous fragment despite x number of timepoints continuously being below the speed threshold 
    vf = processed_behaviourData.vel_for;
    vs = processed_behaviourData.vel_side;
    speed = sqrt(vf.^2 + vs.^2); 
    
    [index] = find(speed >= speedThres); 
    count = 1;
    shiftPoint = [];
    for i = 2:length(index)
        if index(i) - index(i-1) > contThres 
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
     
    no0Vel_idx = index; 
    no0Vel_frag = frag; 

end
