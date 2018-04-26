function out = plot_xtick(self, nticks)
if nargin<2
    nticks = floor(self.para.duration)+1;
end

out = linspace(0, self.para.trial_length, nticks);
if self.include_iti
    iti_length = self.para.iti_length;
    out = [0, out + iti_length, self.sample_length];
end
end