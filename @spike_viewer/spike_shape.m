function [spikeShape, spikeShapeAverage, spike_duration_ms, xplot, single_spike_shape] = spike_shape(self, pre_ms, pos_ms)

if nargin<3
    pos_ms = 0; % extract recording 0 ms after refractory period
end

if nargin<2
    pre_ms = 0; % extract recording 1 ms before start of spike
end

rec3d = self.rec3d_repeat_para_sample;
spike_trace = self.spike_trace;
spike_length = round((self.refractory_period+pre_ms+pos_ms)/1000*self.para.sample_rate);
xplot = linspace(-pre_ms, self.refractory_period+pos_ms, spike_length)/1000*self.para.sample_rate;
single_spike_shape = cell(self.para.num_repeat, self.para.num_unique_stm);

for repeat = 1:self.para.num_repeat
    for i_stim = 1:self.para.num_unique_stm
        fprintf('SpikeViewer: extracting spike shape, repeat: %d, stimuli: %d\n', repeat, i_stim)
        rec = squeeze(rec3d(repeat,i_stim,:));
        spiketrace = squeeze(spike_trace(repeat,i_stim,:));
        spike_shape_ind = find(spiketrace) - pre_ms/1000*self.para.sample_rate;
        num_spike = length(spike_shape_ind);
        spike_shape_ind = bsxfun(@plus, repmat(spike_shape_ind(:), 1, spike_length), 0:spike_length-1);
        max_ss_ind = max(max(spike_shape_ind));
        length_rec = length(rec(:));
        
        if max_ss_ind>length_rec
            warning('SpikeViewer: spike_index (max: %d) is larger than length of recording (%d)\n', max_ss_ind, length_rec)
            fprintf('Spikeviewer: adding nans to the end of recording\n')
            rec2 = [rec; nan(max_ss_ind-length_rec,1)];
            temp_spike_shape = reshape(rec2(spike_shape_ind), num_spike, spike_length);
        else
            temp_spike_shape = reshape(rec(spike_shape_ind), num_spike, spike_length);
        end
        single_spike_shape{repeat, i_stim} = temp_spike_shape;
    end
    
    spikeShape = cell2mat(single_spike_shape(:))';
    
end

if self.spike_peak_alignment
    alignSpikeShape = nan(size(spikeShape));
    num_all_spike = size(spikeShape,2);
    [~, idx] = max(spikeShape);
    min_idx = min(idx);
    range_idx = range(idx);
    
    align_length = spike_length - range_idx;
    align_idx_1st = (idx - min_idx + 1) + (0:spike_length:numel(spikeShape)-spike_length);
    align_idx = bsxfun(@plus, repmat(align_idx_1st, align_length,1), (0:align_length-1)');
    
    alignSpikeShape(1:align_length,:) = reshape(spikeShape(align_idx(:)), align_length, num_all_spike);
    
    spikeShape = alignSpikeShape;
end

spikeShapeAverage = nanmean(spikeShape,2);

rest_start_ind = round(spike_length*2/3);
spike_post_mean = abs(nanmean(spikeShapeAverage(rest_start_ind:end)));
spike_post_std = nanstd(spikeShapeAverage(rest_start_ind:end));

spike_zero_ind = find(abs(spikeShapeAverage - spike_post_mean)<spike_post_std);
spike_plain_ind = find(abs(diff(spikeShapeAverage))<spike_post_std)+1;
spike_end_ind = intersect(spike_plain_ind, spike_zero_ind);
[~, spike_valley_ind] = min(spikeShapeAverage);
spike_end_ind(spike_end_ind<spike_valley_ind) = [];

if isempty(spike_end_ind)
    warning('SpikeViewer: we have a problem computing the spike_duration\n')
    fprintf('SpikeViewer: spike_valley_ind:%d\n', spike_valley_ind)
    fprintf('SpikeViewer: spike_zero_ind:%d\n', spike_zero_ind)
    fprintf('SpikeViewer: spike_plain_ind:%d\n', spike_plain_ind)
    fprintf('SpikeViewer: you may need to check the spike shape in the plot and change spike_thresh\n')
    
    spike_duration_ms = nan;
else
    spike_duration_ms = spike_end_ind(1)/self.para.sample_rate*1000 - pre_ms;
    spikeShape(:,spike_end_ind(1)+1:end) = [];
    spikeShapeAverage(spike_end_ind(1)+1:end) = [];
    xplot(spike_end_ind(1)+1:end) = [];
end

end