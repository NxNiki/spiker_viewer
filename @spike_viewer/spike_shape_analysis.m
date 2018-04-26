function [spike_range, peak_valley_ratio, peak_valley_duration, recoverSlope, valleyFWHM] = spike_shape_analysis(self, slopeTimePoint, nsteps)
if nargin<2
    slopeTimePoint = 0.5; %ms
end

if nargin<3
    nsteps = 1;
end

[~, average_spike_shape] = self.spike_shape(0,slopeTimePoint);
[spike_peak, spike_peak_ind] = max(average_spike_shape);
[spike_valley, spike_valley_ind] = min(average_spike_shape);

spike_range = spike_peak - spike_valley;
peak_valley_ratio = abs(spike_peak/spike_valley);

peak_valley_duration = (spike_valley_ind - spike_peak_ind)/self.para.sample_rate*1000;

slope_point_ind = slopeTimePoint/1000*self.para.sample_rate + spike_valley_ind;
recoverSlope = (average_spike_shape(slope_point_ind+nsteps) - average_spike_shape(slope_point_ind))/nsteps...
    /(1000/self.para.sample_rate); % mV/ms

average_spike_shape(average_spike_shape>0)=0;
% average_spike_shape(isnan(average_spike_shape))=[];

valleyFWHM = fwhm((1:length(average_spike_shape))/self.para.sample_rate*1000, abs(average_spike_shape));