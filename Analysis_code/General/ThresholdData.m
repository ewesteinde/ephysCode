function [Thres_frag, Thres_idx] = ThresholdData(thres, below, measurement2thres, contThres, processed_trialData)

%% extract timepoint indices when the fly is moving
% at some point ask Jenny for her statisitcal version of this to determine
% what is/isn't a start/stop transition 

% speedThres = minimum speed threshold to define when the fly is and isn't moving, only takes into account side & forward velocities
% contThres = threshold for something to count as a continuous fragment despite x number of timepoints continuously being below the speed threshold 
if below == 1
    measurement2thres = -(measurement2thres);
    thres = -(thres); 
end

Thres_frag = cell(1,length(processed_trialData));
Thres_idx = cell(1,length(processed_trialData));

for t = 1:length(processed_trialData)
    [index] = find(measurement2thres >= thres); 
    count = 1;
    frag = {};
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
     
    Thres_idx = index; 
    Thres_frag{t} = frag; 
end

end