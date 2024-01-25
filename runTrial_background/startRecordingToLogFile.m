% start logging data continuously in the background
function startRecordingToLogFile(niIO)
    if niIO.IsRunning
        niIO.stop();
    end
    niIO.startBackground();
end

