function out = group_mean(self, data, para_ind)

group_var = self.para.unique_condition_prop_val;
data = [mean(data,2), group_var];
out = grpstats(data, group_var(:,para_ind), 'mean');
end