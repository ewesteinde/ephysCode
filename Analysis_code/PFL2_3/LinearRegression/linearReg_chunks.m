function [BVm, ave_deg, ave_rho, prefHeads, R2_Vm] = linearReg_chunks(seg, processed_behaviourData, processed_trialData, CL_chunks, fileName)


count = 0; 
for chunk = seg
    count = count + 1; 
    X = [];
    start = CL_chunks(chunk, 1);
    finish = CL_chunks(chunk, 2);
    
    f = fieldnames(processed_behaviourData);
    for k=1:numel(f)
        chunk_behaviourData.(f{k}) = processed_behaviourData.(f{k})(start:finish); 
    end


    f = fieldnames(processed_trialData);
    for k=1:numel(f)
        chunk_trialData.(f{k}) = processed_trialData.(f{k})(start:finish); 
    end

    prefHeads = [-180:10:170];
    for p = 1:length(prefHeads)
        prefHead = prefHeads(p);
        
    y_Vm = chunk_trialData.scaledOutput_down';
        

        X(:,1) = (cosd(chunk_behaviourData.angle-prefHead) - min(cosd(chunk_behaviourData.angle-prefHead)))/(max(cosd(chunk_behaviourData.angle-prefHead))-min(cosd(chunk_behaviourData.angle-prefHead)));
        % implementing mean normalization
        X(:,2) = (chunk_behaviourData.vel_for - min(chunk_behaviourData.vel_for))/(max(chunk_behaviourData.vel_for) - min(chunk_behaviourData.vel_for)); 

        X(:,3) =(chunk_behaviourData.vel_yaw - min(chunk_behaviourData.vel_yaw))/(max(chunk_behaviourData.vel_yaw) - min(chunk_behaviourData.vel_yaw)); 

    % Lasso/elastic net regularization, must calc R^2 independently
        [B_Vm,FitInfo_Vm] = lasso(X, y_Vm, 'NumLambda',5, 'Alpha', 0.5) ;

        % Coefficients for heading & velocities
        BVm_chunk(:,p) = B_Vm(:,1); 
 
        yCalc_Vm = FitInfo_Vm.Intercept(1) +  X*B_Vm(:,1); 
        
        R2_Vm_chunk(p) =  1 - sum((y_Vm - yCalc_Vm).^2)/sum((y_Vm - mean(y_Vm)).^2);
    end
    
    

    
    R2_Vm_posTheta = R2_Vm_chunk(BVm_chunk(1,:) >= 0); 
    
    % set the chunk index to be the one where both the Bthetas are >= 0 and
    % R2 is the maximal value obtained across fitting the tested preferred headings. Controls for situations where two values of
    % preferred heading maximize R2 where they're 180 deg apart
    
    B_idx_Vm = BVm_chunk(1,:) >= 0 & R2_Vm_chunk == max(R2_Vm_posTheta); 

    R2_Vm(count) = R2_Vm_chunk(B_idx_Vm);

    BVm(1,count) = BVm_chunk(1,B_idx_Vm); 
    BVm(2,count) = BVm_chunk(2,B_idx_Vm); 
    BVm(3,count) = BVm_chunk(3,B_idx_Vm); 
    prefHead_Vm(count) = prefHeads(B_idx_Vm);

            
    [no0Vel_frag, no0Vel_idx] = remove0velocity(1.5, 1, chunk_behaviourData);
    %calculate percent time spent moving for each chunk
    percent_moving(count) = length(no0Vel_idx)/length(chunk_behaviourData.vel_for)*100;
    %calculate average euclidean speed for each chunk
    meanSpeed(count) = median(sqrt(chunk_behaviourData.vel_for(no0Vel_idx).^2 + chunk_behaviourData.vel_side(no0Vel_idx).^2));
    meanVf(count) = median(chunk_behaviourData.vel_for(no0Vel_idx)); 
    

%simply gather all timepoints during which the fly's speed was above the
%threshold

       f = fieldnames(chunk_behaviourData);
    for k=1:numel(f)
        no0vel_CL.(f{k}) = chunk_behaviourData.(f{k})(no0Vel_idx); 
    end
    
           f = fieldnames(chunk_trialData);
    for k=1:numel(f)
        no0vel_CL.(f{k}) = chunk_trialData.(f{k})(no0Vel_idx); 
    end
    
% calc heading vector & strength
    window = 60; %seconds
    minVel = 2; %mm/s
    sampRate = 1000; %Hz 

    [rho, theta, idx_windows, i, k] = PFL2_3_headingDist(window, minVel, chunk_behaviourData, sampRate, fileName, 0);
    
    % vector strength of heading
    ave_rho(count) = median(rho); 
    % circular mean of heading
    ave_theta(count) = median(theta); 
    ave_deg = rad2deg(ave_theta); 
    
    prefHeads = wrapTo180(prefHead_Vm);
    
end
end