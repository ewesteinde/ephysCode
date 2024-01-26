function [N,centers, mean_bin] = create_binned_mean(inputs, outputs, edges)
%create_binned_mean Bins the inputs and then takes the mean of the outputs
%in each bin.
    [N, edges, bin] = histcounts(inputs', edges);
    % Add 1 to the bin_yaw values, but remember to remove the first bin
    % after accumarray
    % bin & output arrays paired together as if in a 2 column array. for
    % each bin(i),output(i) pair it sums the output values of all the pairs
    % with the same bin value, take the mean of this (sum Vm data points
    % for each panel to find ave Vm/panel & so estimate angle at which Vm
    % is at the highest mag --> pref angle)
    temp = accumarray(bin+1, outputs', [length(edges) 1]);
    % divides the summed Vm values in temp by the number of Vm samples held
    % in N --> get mean Vm in each bin
    mean_bin = bsxfun(@rdivide, temp(2:end), N');
    centers = edges(1:end-1)+diff(edges)/2;
end