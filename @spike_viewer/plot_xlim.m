function out = plot_xlim(self, option)
if nargin < 2
    option = 1;
end
switch option
    case 1
        out = [0, self.sample_length];
    case 2
        out = [0, self.sample_length] - self.para.iti_length;
end