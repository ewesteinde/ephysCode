% after you've stopped the session (triggered by button app) this will read in the logged data in a
% mat file, figure out where the options input comes from, the button?
function dataLog = loadDaqDataFromLogFile(dataLog)
    % Read in the data from the logfile
    fclose(dataLog.temp_log_file_id);
    dataLog.temp_log_file_id = fopen(dataLog.temp_log_filepath,'r');
    tmp_data = fread(dataLog.temp_log_file_id,'double');
    %log file is read in as one single column of data, need to split into
    %indivdual DAQ channels
    tmp_data = reshape(tmp_data,14,[]);
    dataLog.daq_data = tmp_data(2:end,:)';
    clear tmp_data % profiler suggests bad cleanup on this var
end

