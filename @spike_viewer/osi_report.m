function [osi, para] = osi_report(self, data_type)

% prepare data:
if strcmp(data_type, 'spike')
    HDdata = self.spike_rate_hyper_array;
    dur = self.para.duration;
else
    HDdata = self.membrane_potential_hyper_array;
    dur = 1;
end
DataMean = mean(HDdata,1)/dur;
clear HDdata;

condition_prop_len = self.para.condition_prop_length;
numtf = condition_prop_len(3);
numsf = condition_prop_len(2);

osi = nan(numtf*numsf,1);
para= nan(numtf*numsf,2);

para_ind= combination([numsf, numtf]);
sf_list = self.para.condition_prop_list.spatial_frequency;
tf_list = self.para.condition_prop_list.temporal_frequency;
para(:,1) = sf_list(para_ind(:,1));
para(:,2) = tf_list(para_ind(:,2));

for i=1:numtf*numsf
    i_sf = para_ind(i,1);
    i_tf = para_ind(i,2);
    data_ori = squeeze(SelectDim(DataMean, [self.dim_info.tf, self.dim_info.sf], [i_tf, i_sf]));
    osi(i) = self.compute_osi(self.para.condition_prop_list.direction, data_ori);
end
end