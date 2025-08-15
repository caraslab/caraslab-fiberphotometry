function filtered_data = doublefilter_FPdata(raw_data, sampling_rate, fslow)
    [b1, a1] = butter(3, fslow / sampling_rate * 2, 'low'); 
    
    % Add padding to remove filtering artifact at the end
    n_pad = floor(5*sampling_rate);
    orig_length = length(raw_data);
    buff = [raw_data; repmat(raw_data(orig_length), n_pad, 1)];

    filtered_data = filter(b1, a1, buff); % causal forward filter
    filtered_data = flipud(filtered_data); % reverse time
    filtered_data = filter(b1, a1, filtered_data); % causal forward filter again
    filtered_data = flipud(filtered_data); % reverse time back
    filtered_data = filtered_data(1:orig_length);  % Remove padding zeros
end