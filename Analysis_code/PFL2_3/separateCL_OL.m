function [CL, OL, CL_startStopIdx, OL_startStopIdx] = separateCL_OL(startOL, endOL, lengthCL, processed_behaviourData, processed_trialData)  %% separate idx at which fly is moving in CL vs OL
    
% start = 1819.45;
% finish = 1844.47;
%     startOL  = find(processed_behaviourData.time == start); 
%     endOL = find(processed_behaviourData.time == finish); 

    lengthOL =  ((endOL/1000)-(startOL/1000))*1000; 
    vf = processed_behaviourData.vel_for;
    
    idx = 1; 
    idx_CL = [];
    idx_OL = []; 
    CL_startStopIdx = [];
    OL_startStopIdx = [];
    count = 1; 
    
    while idx <= length(vf)
        if idx == 1
            tempCL = 1:startOL; 
            tempOL = max(tempCL)+1:max(tempCL)+lengthOL;
            CL_startStopIdx(count,:) = [tempCL(1),tempCL(end)];
            OL_startStopIdx(count,:) = [tempOL(1),tempOL(end)];
            count = count + 1; 
        elseif idx + lengthCL + lengthOL > length(vf)
            if idx + lengthCL >= length(vf)
                tempCL = idx:length(vf);
                tempOL = []; 
                CL_startStopIdx(count,:) = [tempCL(1),tempCL(end)];
                count = count + 1; 
            else
                tempCL = idx:idx+lengthCL-1; 
                tempOL = max(tempCL)+1:length(vf);
                CL_startStopIdx(count,:) = [tempCL(1),tempCL(end)];
                OL_startStopIdx(count,:) = [tempOL(1),tempOL(end)];
                count = count + 1; 
            end
        else 
            tempCL = idx:idx+lengthCL-1; 
            tempOL = max(tempCL)+1:max(tempCL)+lengthOL; 
            CL_startStopIdx(count,:) = [tempCL(1),tempCL(end)];
            OL_startStopIdx(count,:) = [tempOL(1),tempOL(end)];  
            count = count + 1; 
        end
        idx = max(max(tempOL), max(tempCL))+1; 
        idx_CL = [idx_CL,tempCL];
        idx_OL = [idx_OL, tempOL];
    end
    idx_CL = round(idx_CL); 
    idx_OL = round(idx_OL);
CL.vel_for = processed_behaviourData.vel_for(round(idx_CL)); 
CL.vel_side = processed_behaviourData.vel_side(idx_CL); 
CL.vel_yaw = processed_behaviourData.vel_yaw(idx_CL); 
CL.angle = processed_behaviourData.angle(idx_CL); 

CL.scaledOutput_down = processed_trialData.scaledOutput_down(idx_CL); 
CL.smooth_Vm = processed_trialData.smooth_Vm(idx_CL); 
CL.fRate_sec = processed_trialData.fRate_sec(idx_CL); 

OL.vel_for = processed_behaviourData.vel_for(idx_OL); 
OL.vel_side = processed_behaviourData.vel_side(idx_OL); 
OL.vel_yaw = processed_behaviourData.vel_yaw(idx_OL); 
OL.angle = processed_behaviourData.angle(idx_OL); 

OL.scaledOutput_down = processed_trialData.scaledOutput_down(idx_OL); 
OL.smooth_Vm = processed_trialData.smooth_Vm(idx_OL); 
OL.fRate_sec = processed_trialData.fRate_sec(idx_OL); 
end