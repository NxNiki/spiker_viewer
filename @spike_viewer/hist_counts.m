function [counts, edges] = hist_counts(self)
% do hiscounts for each trial:
spike_trace = self.spike_trace;
numRepeat = self.para.num_repeat;
numUniqueStm = self.para.num_unique_stm;
numBins = self.num_bins;

counts = zeros(numRepeat, numUniqueStm, numBins);
edges = linspace(0, self.para.duration, self.num_bins+1);
for i = 1:numRepeat
    for j = 1:numUniqueStm
        try
            counts(i,j,:) = histcounts(find(spike_trace(i,j,:))/self.para.sample_rate, edges);
        catch
            counts(i,j,:) = histc(find(spike_trace(i,j,:))/self.para.sample_rate, edges(1:end-1));
        end
    end
end
end