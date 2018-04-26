function out = plot_xticklabel(self, nticks)
if nargin<2
    nticks = floor(self.para.duration)+1;
end

out = linspace(0, self.para.duration, nticks);
if self.include_iti
    out = [-self.para.iti, out, self.sample_duration-self.para.iti];
end
end