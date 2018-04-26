function out = plot_xtime(self)

if self.include_iti
    out = linspace(-self.para.iti,self.para.duration+self.para.iti, self.sample_length);
else
    out = (1:self.para.duration)/self.para.sample_rate;
end

end