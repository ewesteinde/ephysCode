function [chunks] = findContinuousData(data)

    count = 1; 
    prev_jump = 0; 
    transition = zeros(size(data)); 

    for idx = 1:length(data)
        value = data(idx);
        if idx == 1 
            chunks(count,1) = idx;
            prev_val = value; 
        elseif idx == length(data)
            chunks(count,2) = idx; 
        else
            if value ~= prev_val && ~(isnan(value)) && value == 0 
                chunks(count,1) = idx;
            elseif value ~= prev_val && ~(isnan(value)) && value == 1 
                chunks(count,2) = idx;
                count = count + 1; 
            end
            prev_val = value; 
        end
    end
end