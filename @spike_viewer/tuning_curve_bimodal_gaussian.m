function f = tuning_curve_bimodal_gaussian(self, data_type, manual_pref)
% bimodal gaussian fit of tuning curve. data_type is 'spike' or 'mp'.
% manual_pref is a n by 3 matrix. each row corresponds to a combination of
% temporal frequency and spatial frequency that you want to manual set the
% preferred orientaton. 1st column is the temporal frequency idx, 2nd colum
% is the spatial frequency idx, 3rd is the direction you want to manually
% set.
% by Niki 2016/4/7.


% make sure the condition property is set with the right order: direction,
% spatial frequency, temporal freqency
condition_prop_len = self.para.condition_prop_length;
numtf = condition_prop_len(3);
numsf = condition_prop_len(2);
numdir= condition_prop_len(1);

manual_pref_idx = nan(numtf, numsf);
if nargin==3
    ind = sub2ind([numtf, numsf], manual_pref(1,:), manual_pref(2,:));
    manual_pref_idx(ind) = manual_pref(3,:);
end

if nargin<2
    data_type='spike';
end

self.spike_rate_minus_spontaneuous = true;

% prepare data:
if strcmp(data_type, 'spike')
    HighDimData = self.spike_rate_hyper_array;
    dur = self.para.duration;
else
    HighDimData = self.membrane_potential_hyper_array;
    dur = 1;
end
DataMean = mean(HighDimData,1)/dur;
clear HighDimData;

% fit function:
bigaussian = @(R0,A1,A2,preftheta,sigma,theta)(...
    R0 + ...
    A1*exp(1).^((cos(theta-preftheta)-1)/sigma^2) + ...
    A2*exp(1).^((cos(theta-preftheta-pi)-1)/sigma^2)...
    );

fittobj = fittype(bigaussian,'coefficients',{'R0','A1','A2','preftheta','sigma'},'independent','theta');

% prepare parameter information:
tf_list = self.para.condition_prop_list.temporal_frequency;
sf_list = self.para.condition_prop_list.spatial_frequency;
dir_list = self.para.condition_prop_list.direction;



% plot:
f = figure('position', [0, 0, numsf*400, numtf*300]);
lim = [min(DataMean(:))-0.1, max(DataMean(:))+0.1];
theta = 0:.1:2*pi;

for i_tf = 1: numtf
    for i_sf = 1:numsf
        plot_ind = (i_tf - 1) * numsf + i_sf;
        data_ori = squeeze(SelectDim(DataMean, [self.dim_info.tf, self.dim_info.sf], [i_tf, i_sf]));
        OSI = self.compute_osi(self.para.condition_prop_list.direction, data_ori);
        axe = subplot(numtf, numsf, plot_ind);
        
        
        if max(dir_list)>2*pi
            degree2pi = 180/pi;
        else
            degree2pi = 1;
        end
        
        R0 = self.spontaneous_spike_rate;
        [A1 ,pref_idx] = max(data_ori(:));
        sigma = std(data_ori(:));

        
        fit_model = cell(1,length(pref_idx));
        gof_table = nan(1,length(pref_idx));
        for i_pref_dir = 1:length(pref_idx);
            i_dir = pref_idx(i_pref_dir);
            opposite_dir_idx = mod(i_dir + numdir/2, numdir);
            A2 = data_ori(opposite_dir_idx);
            preftheta = dir_list(i_dir)/degree2pi;
            % least square fit of bigaussian:
            [fit1,gof,fitinfo] = fit(dir_list/degree2pi,data_ori(:),fittobj,'StartPoint',[R0,A1,A2,preftheta,sigma]);
            disp(fitinfo.message)
            fit_model{i_pref_dir} = fit1;
            gof_table(i_pref_dir) = gof.sse;
        end
        
        [~,minsse_idx] = min(gof_table);
        fit1 = fit_model{minsse_idx};
        fitcurve = fit1(theta);
        plot(axe, theta, fitcurve, 'b--', 'LineWidth', 1, 'MarkerSize', 1);
        hold on
        plot(axe, dir_list/degree2pi, data_ori, '*', 'LineWidth', 1, 'MarkerSize', 2);
        hold off
        ylim(lim)
        
        title(sprintf([data_type,'OSI = %0.3f'], OSI));
        set(gca,'box','off','XTick',self.para.condition_prop_list.direction);
        
        if i_tf == numtf
            xlabel(axe, sprintf('spatial frequency: %0.2f', sf_list(i_sf)),'FontSize',12)
        end
        if i_sf == 1
            ylabel(axe, sprintf('temporal frequency: %0.2f', tf_list(i_tf)),'FontSize',12)
        end
    end
end
end