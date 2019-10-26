function [bin_data] = extract_data(raw_data,f)
%EXTRACT_DATA Extract binned data from timespike values
%   raw_data - raw timespike data
%   f - bin frequency 
mid_data = [];
for n = 1:length(raw_data.S.C)
    mid_data{n} = raw_data.S.C{n}.tsd.t;
end
end_data.spike_times = [];
for idx = 1:length(mid_data)
    end_data(idx).spike_times = mid_data{idx};
end
min = end_data(1).spike_times(1);
for j = 1:length(end_data(1).spike_times)
    if end_data(1).spike_times(j) < min
        min = end_data(1).spike_times(j);
    end
end
max = end_data(1).spike_times(1);
for j = 1:length(end_data(1).spike_times)
    if end_data(1).spike_times(j) > max
        max = end_data(1).spike_times(j);
    end
end
bin_data = binspikes(end_data,f,[min,max]);
end
