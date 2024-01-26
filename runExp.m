if ~exist('date','var')
    ephysSettings
    rootPath = settings.mainDataDir;
    date = input('Date: ' ,'s');
    foldername = fullfile(rootPath,'newData',date);
    cell = 0;
    nfly = 0;
end
%% Call every new trial

clearvars t_*
newCell = input('New cell? y/n ','s');
if newCell == 'y'
    cell = cell + 1;
    trial= 1;
    newFly = input('New fly? y/n ','s');
    if newFly == 'y'
        nfly = nfly + 1;
        fly.flyExp = input('Brief experiment descrption: ', 's');
        newGenotype = input('new genotype? y/n ','s');
        if newGenotype == 'y'
             fly.flyGenotype = input('Fly Genotype? ', 's');
        end
        fly.flyAge = input('Fly Age? ', 's');
        fly.flySex = input('Fly Sex? ', 's');
        fly.flyNotes = input('Fly Notes: ', 's');
        fly.timestamp = datetime;
        fly.flyName = strcat(date,'_',fly.flyGenotype,'_0',num2str(nfly));
    end 
end

trialnum = strcat('trial_',int2str(trial));
cellnum = strcat('cell_',int2str(cell));

% Choose acquisition type and collect trace

trialType = input('sim, acc, pip, iR, cPulse, lPulse, dark, jump, CL, OL, CL_OL, bjumpCL, jumpCL, ionto_pCL, ionto_CL, ionto_sCL, ionto : ','s');

if strcmp(trialType, 'pip')|| strcmp(trialType,'acc')
    check = input('Seal test on & Vclamp? y/n ','s');
    if strcmp(check, 'y')
        if strcmp(trialType, 'pip')
            [pipetteR, PipetteRData, PipetteRMeta, PipetteRDataraw] = measurePipetteResistance();
        elseif strcmp(trialType,'acc')
            seal = input('seal? y/n ','s');
            if strcmp(seal, 'y')
                [holdingCurrent, accessResistance, sealResistance, SealRData, SealRMeta, SealRDataraw] = measureAccessResistance();
            else
                [holdingCurrent, accessResistance, inputResistance, AccessRData, AccessRMeta, AccessRDataraw] = measureAccessResistance();
            end
        elseif strcmp(trialType,'iR')
            check = input('Ext input on & Iclamp? y/n ','s');
            if strcmp(check, 'y')
                [t_inputData,t_trialInR] = measureInputResistance(); % TODO: update for JL rig
            end
        end
    end
else
    if strcmp(trialType, 'sim')
        newParam = input('new parameters? y/n ','s');
        if newParam == 'y'
            trialRepeats = input('Number of trials ');
            trialLength = input('Length of trial (s) ');
            inR = input('inR? 0/1 ');
        end
        if inR == 1
            check = input('Ext input on & Iclamp? y/n ','s');
            if strcmp(check, 'y')   
                [rawData, trialData, trialMeta] = acquireSimpleTrial(trialRepeats, trialLength,inR);
            end
        else
            [rawData, trialData, trialMeta] = acquireSimpleTrial(trialRepeats, trialLength,inR);
        end
    elseif strcmp(trialType, 'cPulse')
        newParam = input('new parameters? y/n ','s');
        if newParam == 'y'
            trialRepeats = input('Number of trials ');
            amp = input('amplitude of stimulus (pA)');
            pulseLength = input('duration of stim (s)');
            bufferLength = input('time between stimuli (s)');
            inR = input('inR? 0/1 ');
        end 
        check = input('Ext input on & Iclamp? y/n ','s');
        if strcmp(check, 'y')
            [rawData, trialData, trialMeta] = acquireCurrentPulseTrial(trialRepeats, amp, pulseLength, bufferLength,inR); 
        end
    elseif strcmp(trialType,'lPulse')
        newParam = input('new parameters? y/n ','s');
        if newParam == 'y'
            trialRepeats = input('Number of trials ');
            pulseLength = input('duration of stim (ms)');
            bufferLength = input('time between stimuli (s)');
            inR = input('inR? 0/1 ');
        end
        filter_check = input('Correct filter? y/n ','s');
        if strcmp(filter_check, 'y')
            [rawData, trialData, trialMeta] = acquireLightPulseTrial(trialRepeats, pulseLength, bufferLength, inR);
        end
    elseif strcmp(trialType, 'CL_OL')
        newParam = input('new parameters? y/n ','s');
        if newParam == 'y'
            trialLength = input('Length of trial (s) ');
        end
            [rawData, trialData, trialMeta, behaviorData, tempFicTrac] = acquireCL_OL_Trial(trialLength);
            trialMeta.pattern = settings.panels.barPattern;  
    elseif strcmp(trialType, 'CL')
        newParam = input('new parameters? y/n ','s');
        if newParam == 'y'
            trialRepeats = input('Number of trials ');
            trialLength = input('Length of trial (s) ');
            inR = 0; %input('inR? 0/1 ');
        end
        if inR == 1
            check = input('Ext input on & Iclamp? y/n ','s');
            if strcmp(check, 'y')   
                [rawData, trialData, trialMeta, tempFicTrac] = acquireCLTrial(trialRepeats, trialLength,inR);
                trialMeta.pattern = settings.panels.barPattern;
            end
        else
            [rawData, trialData, trialMeta, behaviorData, tempFicTrac] = acquireCLTrial(trialRepeats, trialLength, inR);
            trialMeta.pattern = settings.panels.barPattern;
        end
    elseif strcmp(trialType, 'dark')
        newParam = input('new parameters? y/n ','s');
        if newParam == 'y'
            trialRepeats = input('Number of trials ');
            trialLength = input('Length of trial (s) ');
            inR = input('inR? 0/1 ');
        end
        if inR == 1
            check = input('Ext input on & Iclamp? y/n ','s');
            if strcmp(check, 'y')   
                [rawData, trialData, trialMeta, behaviorData] = acquireDarkTrial(trialRepeats, trialLength,inR);
                trialMeta.pattern = settings.panels.darkPattern;
            end
        else
            [rawData, trialData, trialMeta, behaviorData] = acquireDarkTrial(trialRepeats, trialLength, inR);
            trialMeta.pattern = settings.panels.darkPattern;
        end
    elseif strcmp(trialType, 'jump')
    newParam = input('new parameters? y/n ','s');
    if newParam == 'y'
        trialRepeats = input('Number of trials ');
        trialLength = input('Length of trial (s), jump every 2 min ');
        inR = input('inR? 0/1 ');
    end
    if inR == 1
        check = input('Ext input on & Iclamp? y/n ','s');
        if strcmp(check, 'y')   
            [rawData,trialData, trialMeta, behaviorData] = acquireJumpTrial(trialRepeats, trialLength, inR);
            trialMeta.pattern = settings.panels.jumpPattern;
        end
    else
        [rawData,trialData, trialMeta, behaviorData] = acquireJumpTrial(trialRepeats, trialLength, inR);
        trialMeta.pattern = settings.panels.jumpPattern;
    end
elseif strcmp(trialType, 'OL')
    newParam = input('new parameters? y/n ','s');
    if newParam == 'y'
        trialRepeats = input('Number of trials ');
        trialLength = input('Length of trial (s) ');
        inR = input('inR? 0/1 ');
    end
        
        function_number = input('function number? 19 = jump, 9 = spin R, 7 = spin L, 8 = bar_prefHead ');
            if ~(function_number == 8)
                if function_number == 19
                    pattern_number = 3;
                else
                    pattern_number = 1;
                end
            else
                pattern_number = 5;
            end
        starting_position = input('starting pattern position x? ');
        if inR == 1
            check = input('Ext input on & Iclamp? y/n ','s');
            if strcmp(check, 'y')   
                [rawData,trialData, trialMeta, behaviorData] = acquireOLTrial(trialRepeats, trialLength,inR, pattern_number, function_number, starting_position);
            end
        else
            [rawData, trialData, trialMeta, behaviorData] = acquireOLTrial(trialRepeats, trialLength, inR, pattern_number, function_number, starting_position);
        end
elseif strcmp(trialType,'ionto_pCL')
   newParam = input('new parameters? y/n ','s');
    if newParam == 'y'
        trialRepeats = input('Number of trials ');
        stimLength = input('Ionto pulse length (ms) ');
        bufferLength = input('Time between pulses (s) ');
    end
    [rawData, trialData, trialMeta, behaviorData,tempFicTrac] = acquireIonto_CLPulse_Trial(trialRepeats, stimLength, bufferLength);
elseif strcmp(trialType,'ionto_CL')
   newParam = input('new parameters? y/n ','s');
    if newParam == 'y'
        trialLength = input('Length of trial (s) ');
        stimLength = input('Ionto pulse length (ms) ');
        bufferLength = input('Time between pulses (s) ');
    end
    [rawData, trialData, trialMeta, behaviorData,tempFicTrac] = acquireIonto_CL_Trial(trialLength, stimLength, bufferLength);
elseif strcmp(trialType,'ionto')
   newParam = input('new parameters? y/n ','s');
    if newParam == 'y'
        trialRepeats = input('Number of trials ');
        stimLength = input('Ionto pulse length (ms) ');
        bufferLength = input('Time between pulses (s) ');
    end
    [rawData, trialData, trialMeta, behaviorData] = acquireIontoPulseTrial(trialRepeats, stimLength, bufferLength);
elseif strcmp(trialType,'ionto_sCL')
       newParam = input('new parameters? y/n ','s');
    if newParam == 'y'
        trialLength = input('Length of trial (s) ');
        bufferLength = input('Time between pulses (s) ');
    end
    [rawData, trialData, trialMeta, behaviorData, tempFicTrac] = acquireIonto_CLsetStim_Trial(trialLength, bufferLength);
elseif strcmp(trialType,'jumpCL')
       newParam = input('new parameters? y/n ','s');
    if newParam == 'y'
        trialLength = input('Length of trial (s) ');
    end
    [rawData, trialData, trialMeta, behaviorData, tempFicTrac] = acquire_CL_jump_Trial(trialLength);
 elseif strcmp(trialType,'bjumpCL')
       newParam = input('new parameters? y/n ','s');
    if newParam == 'y'
        trialLength = input('Length of trial (s) ');
        trialRepeats = input('Number of trials ');
    end   
    [rawData, trialData, trialMeta, behaviorData, tempFicTrac] = acquire_CLbehaviour_jump_Trial(trialRepeats, trialLength);
end

end


% Choose to save trial or not
pause(1); 
if ~strcmp(string(trialType), 'pip') && ~strcmp(string(trialType), 'acc') && ~strcmp(trialType, 'iR')
    toSave = input('save? y/n ','s');
    if strcmp(toSave,'y')
        
        filename = fullfile(foldername,['fly_',num2str(nfly)],cellnum,[trialType,'_',trialnum]); 
         
        if ~exist(filename, 'dir')
            mkdir(filename)
        end
        
        if exist('tempFicTrac','dir')
            copyfile(tempFicTrac,filename)
            rmdir(tempFicTrac,'s')
        end             
    
        notes = input('add note? y/n ','s');
            if strcmp(notes, 'y')   
                trialMeta.notes = input('write note: ','s');
            end
        if ~strcmp(trialType,'bjumpCL')    
            finalAcc = input('measure final AccessR? y/n ','s');
            if strcmp(finalAcc, 'y')   
                [holdingCurrentEnd, accessResistanceEnd, inputResistanceEnd, AccessRDataEnd, AccessRMetaEnd, AccessRDatarawEnd] = measureAccessResistance();
                trialMeta.holdingCurrentEnd = holdingCurrentEnd;
                trialMeta.accessREnd = accessResistanceEnd;
                trialMeta.inputREnd = inputResistanceEnd; 
            end
        end
        
        dataPath = fullfile(filename,'trialData.mat');
        rawPath = fullfile(filename, 'rawData.mat');
        metaPath = fullfile(filename,'trialMeta.mat');
        if ~strcmp(trialType,'bjumpCL')
            trialMeta.sealR = sealResistance;
        end
        trialMeta.trialType = trialType;
        
        if strcmp(trialMeta.trialType,'OL')
            if function_number == 8
                trialMeta.funcType = 'heading fixed';
            elseif function_number == 7
                trialMeta.funcType = 'spin left';
            elseif function_number == 9
                trialMeta.funcType = 'spin right';
            end
        end
        
        if ~strcmp(trialType,'bjumpCL')
            trialMeta.pipetteR = pipetteR;
            trialMeta.holdingCurrent = holdingCurrent;
            if exist('accessResistance')
                trialMeta.accessR = accessResistance;
            end
            if exist('inputResistance')
                trialMeta.inputR = inputResistance;  %need to deal with inR from simple vs from access
            end
        end
        trialMeta.settings = settings;
        trialMeta.fly = fly;
        
        %EphysExpSum(trialMeta, date, cell, trial) %adds meta data to
        %onedrive excel sheet for easy reference - 2022/07/18 file got
        %corrupted, can't open anymore
        
        save(dataPath, 'trialData','-v7.3');
        save(metaPath, 'trialMeta','-v7.3');
        save(rawPath, 'rawData','-v7.3'); % -v7.3 creates files slightly larger than the default save option but 1.6x as fast. If you need to go even faster use the '-nocompression' option, files will be very large though
        
        
        if ~strcmp(trialType,'bjumpCL')
            dataPath = fullfile(filename,'PipetteRData.mat');
            rawPath = fullfile(filename, 'PipetteRDataraw.mat');
            metaPath = fullfile(filename,'PipetteRMeta.mat');
            PipetteRMeta.trialType = trialType;
            PipetteRMeta.pipetteR = pipetteR;
            PipetteRMeta.settings = settings;
            PipetteRMeta.fly = fly;
        

            save(dataPath, 'PipetteRData','-v7.3');
            save(metaPath, 'PipetteRMeta','-v7.3');
            save(rawPath, 'PipetteRDataraw','-v7.3');


            dataPath = fullfile(filename,'SealRData.mat');
            rawPath = fullfile(filename, 'SealRDataraw.mat');
            metaPath = fullfile(filename,'SealRMeta.mat');
            SealRMeta.sealR = sealResistance;
            SealRMeta.trialType = trialType;
            SealRMeta.pipetteR = pipetteR;
            SealRMeta.settings = settings;
            SealRMeta.fly = fly;

            save(dataPath, 'SealRData','-v7.3');
            save(metaPath, 'SealRMeta','-v7.3');
            save(rawPath, 'SealRDataraw','-v7.3');


            dataPath = fullfile(filename,'AccessRData.mat');
            rawPath = fullfile(filename, 'AccessRDataraw.mat');
            metaPath = fullfile(filename,'AccessRMeta.mat');
            AccessRMeta.sealR = sealResistance;
            AccessRMeta.trialType = trialType;
            AccessRMeta.pipetteR = pipetteR;
            AccessRMeta.holdingCurrent = holdingCurrent;
            AccessRMeta.accessR = accessResistance;
            if exist('inputResistance')
                AccessRMeta.inputR = inputResistance;  %need to deal with inR from simple vs from access
            end
            AccessRMeta.settings = settings;
            AccessRMeta.fly = fly;

            save(dataPath, 'AccessRData','-v7.3');
            save(metaPath, 'AccessRMeta','-v7.3');
            save(rawPath, 'AccessRDataraw','-v7.3');

            if strcmp(finalAcc, 'y')   
                dataPath = fullfile(filename,'AccessRDataEnd.mat');
                rawPath = fullfile(filename, 'AccessRDatarawEnd.mat');
                metaPath = fullfile(filename,'AccessRMetaEnd.mat');
                AccessRMetaEnd.sealR = sealResistance;
                AccessRMetaEnd.trialType = trialType;
                AccessRMetaEnd.pipetteR = pipetteR;
                AccessRMetaEnd.holdingCurrent = holdingCurrentEnd;
                AccessRMetaEnd.accessR = accessResistanceEnd;
                AccessRMetaEnd.inputR = inputResistanceEnd;  %need to deal with inR from simple vs from access
                AccessRMetaEnd.settings = settings;
                AccessRMetaEnd.fly = fly;

                save(dataPath, 'AccessRDataEnd','-v7.3');
                save(metaPath, 'AccessRMetaEnd','-v7.3');
                save(rawPath, 'AccessRDatarawEnd','-v7.3');
            end
        end
        
        if ~strcmp(string(trialType), 'sim')
            behaviorPath = fullfile(filename,'behaviorData.mat');
            save(behaviorPath, 'behaviorData','-v7.3');
        end
        
        if strcmp(string(trialType), 'CL') || strcmp(string(trialType), 'dark') 
            behaviorPath = fullfile(filename,'behaviorData.mat');
            %if matlab says file too big add '-v7.3' arg to end of save
            save(behaviorPath, 'behaviorData','-v7.3');
        end
        
        trial = trial + 1;
    end    
end
    
        
        
 
