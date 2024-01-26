function logDaqData(~,evt,fid)
    % keep only the inds specified (if specified)
    data = [evt.TimeStamps, evt.Data]';
    % precision is only single
    fwrite(fid, data,'double');
end