function adjustedOffset = adjustOffsetBasedGain(input, gain)
% unusual offsets based on gain settings
    if gain == 0.5
        adjustedOffset = input + 0; 
    elseif gain == 1
        adjustedOffset = input + 0;
    elseif gain == 2
        adjustedOffset = input + 0; 
    elseif gain == 5
        adjustedOffset = input + 0;
    elseif gain == 10
        adjustedOffset = input + 0; 
    elseif gain == 100
        adjustedOffset = input + 0; 
    else
        adjustedOffset = input + 0; 
    end
end