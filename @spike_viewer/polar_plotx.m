function out = polar_plotx(self, prop_name)
list = unique(self.para.condition_prop.(prop_name));
out = [list(:)', list(1)];
end