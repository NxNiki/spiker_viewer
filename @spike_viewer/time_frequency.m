function [f, y, isSimple, para_val] = time_frequency(self, para_condition_prop)

para_val = self.para.unique_condition_prop_val;
if nargin<2
    ind = 1: self.para.num_unique_stm;
else
    ind = self.para.find_condition_prop(para_condition_prop,'unique');
    para_val = para_val(ind,:);
end

numBins = self.num_bins;
psth_data = self.hist_counts;
psth_stml = psth_data(:,ind,:);

frame_per_bin = self.sample_length/numBins;
Fs = self.para.sample_rate/frame_per_bin;                   % Sampling frequency
% T = 1/Fs;                     % Sample time
L = numBins;                   % Length of signal
%t = (0:L-1)*T;                % Time vector
NFFT = L;

Y = fft(squeeze(mean(psth_stml)),NFFT,2)/L;
y = 2*abs(Y(:,1:NFFT/2));

f = Fs/2*linspace(0,1,NFFT/2);

y(:,1) = y(:,1)/2;

F0 = y(:,1);
F1_pos = abs(f-2) == min(abs(f-2));
F1 = y(:,F1_pos);

isSimple = F1./F0;
end