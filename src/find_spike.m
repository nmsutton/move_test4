function spike_found = find_spike(ni, t, firings) 
	% report spike times in a bin of time
	spike_found = false;
	all_spike_times = (find(firings(:,2)==ni));
	for si=1:size(all_spike_times)
		spike_time = firings(all_spike_times(si),1);
        if spike_time == t
            spike_found = true;
        end
    end
end