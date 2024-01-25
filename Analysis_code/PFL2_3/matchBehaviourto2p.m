%% Only for integrating G3 information with fictrac information
function [ftData] = matchBehaviourto2p(folder)
    
%% Get fictrac folder
    
    if strcmp(folder(end),'.')
        folder = folder(1:end-2); 
    end    
    
    ftDir = folder;
    daqFile_info = dir(fullfile(folder,'*_daqData_*.mat'));
    
        %% Get trial data
        expID = get_expID(folder);
        expList = {expID};
        %[expMetadata, trialMetadata, patternMetadata, fictracMetadata] = load_metadata(expList, folder);
        [expMetadata,trialMetadata, patternMetadata, fictracMetadata] = load_metadata(expList, folder);
    
    %% Get variables
    maxVal = 10;
    minVal = 0;
        arenaExtent = 360.;
        initialAngle = -9.375;
        ball = 9;
        luminance = 1;
        yDimxDim = 1;
        cuePosAngleRel = 1;

                
    %% Announce
    disp(['Looking for fictrac data in: ', ftDir])
    
    % Create table of final analysis data
    ftData = [];
    for iTrial = 1:length(daqFile_info)
        load(fullfile(folder,daqFile_info(iTrial).name),'trialData')

        newRow = table({expID}, iTrial, 'VariableNames', {'expID', 'trialNum'});
            
        % Get rates
        fictrac_rate = fictracMetadata.fictracRate; 
        
        try
            NiDaq_rate = trialMetadata.daqSampRate(trialMetadata.trialNum == iTrial);
        catch
            NiDaq_rate = expMetadata.daqSampRate(1); 
        end
        
        maxFlyVelocity = 20; % radians/sec

        % Copy important variables, converting units as needed
        trialData = timetable2table(trialData); 
        
        %% Elena: all Berg4 trials prior to 05/25/22 have the channels for IntSide & IntFor switched
        [ velYaw , intYaw ] = ficTracSignalDecoding_2p( trialData.ficTracYaw , NiDaq_rate , fictrac_rate/2, fictrac_rate, maxFlyVelocity);
        [ velSide , intSide ] = ficTracSignalDecoding_2p( trialData.ficTracIntSide , NiDaq_rate , fictrac_rate/2,fictrac_rate, maxFlyVelocity);
        [ velForward , intForward ] = ficTracSignalDecoding_2p( trialData.ficTracIntForward , NiDaq_rate , fictrac_rate/2, fictrac_rate, maxFlyVelocity);
        
        newRow.trialTime = {seconds(resample_with_padding(seconds(trialData.Time),fictrac_rate,NiDaq_rate))'};        % seconds
        newRow.intFor = {(intForward * ball/2)'};      % mm
        newRow.intSide = {(intSide * ball/2)'};        % mm
        newRow.intHD = {(intYaw)'};                    % radians  
        newRow.velFor = {(velForward * ball/2)'};      % mm/sec
        newRow.velSide = {(velSide * ball/2)'};        % mm/sec
        newRow.velYaw = {(velYaw)'};                   % radians/sec

        % Calculate visual cue position based on X & Y channel pos 
        if ismember('PanelsXDimTelegraph',trialData.Properties.VariableNames) && ismember('PanelsYDimTelegraph',trialData.Properties.VariableNames)
            xframes = patternMetadata.x_num;
            yframes = patternMetadata.y_num;
            pixelAngle = arenaExtent/xframes;

            % calculates width of most salient visual cue, ignores fainter
            % background patterns if present but currently requires main cue to be
            % an individual shape of uniform width & either the brightest or darkest component of the pattern

            pattern_2D = patternMetadata.Pats(:,:,1,1);

            if luminance == 0 
                pattern_1D = pattern_2D(1,:) == min(pattern_2D(1,:));
            else
                pattern_1D = pattern_2D(1,:) == max(pattern_2D(1,:));
            end

            cueWidth = sum(pattern_1D);

            XvoltsPerStep = (maxVal-minVal)./(xframes);
            YvoltsPerStep = (maxVal-minVal)./(yframes);

            % Set limits on voltage
            rawPanelsData = [resample_with_padding(trialData.PanelsXDimTelegraph, fictrac_rate, NiDaq_rate); resample_with_padding(trialData.PanelsYDimTelegraph, fictrac_rate, NiDaq_rate)]';
            rawPanelsData(rawPanelsData < minVal) = minVal;
            rawPanelsData(rawPanelsData > maxVal) = maxVal;

            % Calculate the frame number (round to nearest integer), & calculate the
            % pixel angle of the bar given the frame number.

            frX = round((rawPanelsData(:,1) - minVal)./XvoltsPerStep);
            frY = round((rawPanelsData(:,2) - minVal)./YvoltsPerStep);
            yframeStep = xframes/yframes; 

            % takes into account a change in cue pos due to a yframe
            % step
            % may need to change sign depending on the relationship b/w
            % your x & y frame pos values
            disp('Temporary message for Elena: if processing experiments acquired before 12/16/21 flip yframeStepsign to -');   
            if yDimxDim == 0
                cuePos = mod(frX - yframeStep*(frY-1),xframes);
            else
                cuePos = mod(frX + yframeStep*(frY-1),xframes);
            end
            cuePos(cuePos==0) = xframes;

            % takes into account cue width & starting pos angle
            % relative to fly 
            % may need to change sign depending on the relationship b/w
            % cus pos and and cue angle rel to the fly

            if cuePosAngleRel == 0
                cueAngle = (initialAngle - ((cuePos - 2) + cueWidth/2).*pixelAngle);
            else
                cueAngle = (initialAngle + ((cuePos - 2) + cueWidth/2).*pixelAngle);
            end

            cueAngle = wrapTo180(cueAngle);

            % clean up/filter (artifacts often present at 180 to -180 transitions)   
            cueAngle = smoothdata(cueAngle,'movmedian',5); % can play around with this step
            cueAngle(cueAngle > 180) = 180;
            cueAngle(cueAngle < -180) = -180; 

            % debug plot 
            % figure();plot(cueAngle)

            % Berg4 default: cue pos counterclockwise to the fly = - angles
            %                cue pos clockwise to the fly = + angles 
            % may vary by arena depending on how it/fictrac was setup 
            % IMPORTANT: check your arena's coordinate frame

            newRow.PanelsX = {frX};
            newRow.PanelsY = {frY};
            newRow.cuePos = {cuePos};
            newRow.cueAngle = {cueAngle};

        elseif ismember('PanelsXDimTelegraph',trialData.Properties.VariableNames)
            xframes = patternMetadata.x_num;           
            pixelAngle = arenaExtent/xframes;

            % calculates width of most salient visual cue, ignores fainter
            % background patterns if present but currently requires main cue to be
            % an individual shape of uniform width & brightest component of the pattern

            pattern_2D = patternMetadata.Pats(:,:,1,1);
            pattern_1D = pattern_2D(1,:) == max(pattern_2D(1,:));
            cueWidth = sum(pattern_1D);
            XvoltsPerStep = (maxVal-minVal)./(xframes);


            % Set limits on voltage
            rawPanelsData = [resample_with_padding(trialData.PanelsXDimTelegraph,fictrac_rate, NiDaq_rate)];
            rawPanelsData(rawPanelsData < minVal) = minVal;
            rawPanelsData(rawPanelsData > maxVal) = maxVal;

            % Calculate the frame number (round to nearest integer), & calculate the
            % pixel angle of the bar given the frame number.

            frX = round((rawPanelsData - minVal)./XvoltsPerStep);
            cuePos = frX;

            % takes into account cue width & starting pos angle
            % relative to fly 

            if cuePosAngleRel == 0
                cueAngle = (initialAngle - ((cuePos - 2) + cueWidth/2).*pixelAngle);
            else
                cueAngle = (initialAngle + ((cuePos - 2) + cueWidth/2).*pixelAngle);
            end

            cueAngle = wrapTo180(cueAngle);

            % clean up/filter (artifacts often present at 180 to -180 transitions)

            cueAngle = smoothdata(cueAngle,'movmedian',5); % can play around with this step
            cueAngle(cueAngle > 180) = 180;
            cueAngle(cueAngle < -180) = -180; 

            % debug plot 
            % figure();plot(cueAngle)

            % Berg4 default: cue pos counterclockwise to the fly = - angles
            %                cue pos clockwise to the fly = + angles 
            % may vary by arena depending on how it/fictrac was setup 
            % IMPORTANT: check your arena's coordinate frame

            newRow.PanelsX = frX;
            newRow.cuePos = cuePos;
            newRow.cueAngle = cueAngle;
        end

        % Append to main table
        ftData = [ftData; newRow];
       
    end        
end