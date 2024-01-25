% clears your temporary log file 
function clearLogFileData(dataLog)
    % Close and re-open so that it overwrites content ('w+')
    [~] = fclose(dataLog.temp_log_file_id);
    dataLog.temp_log_file_id = fopen(dataLog.temp_log_filepath,'w+');
end