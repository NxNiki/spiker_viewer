function f = tuning_curve(self, data_type, figure_type)

if nargin<3
    figure_type='normal';
end

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

% prepare parameter information:
tf_list = self.para.condition_prop_list.temporal_frequency;
sf_list = self.para.condition_prop_list.spatial_frequency;

condition_prop_len = self.para.condition_prop_length;
numtf = condition_prop_len(3);
numsf = condition_prop_len(2);

% plot:
f = figure('position', [0, 0, numsf*400, numtf*300]);
lim = [min(DataMean(:))-0.1, max(DataMean(:))+0.1];

for i_tf = 1: numtf
    for i_sf = 1:numsf
        plot_ind = (i_tf - 1) * numsf + i_sf;
        data_ori = squeeze(SelectDim(DataMean, [self.dim_info.tf, self.dim_info.sf], [i_tf, i_sf]));
        OSI = self.compute_osi(self.para.condition_prop_list.direction, data_ori);
        axe = subplot(numtf, numsf, plot_ind);
        
        if strcmp(figure_type,'polar')
            plotx_dir = self.polar_plotx('direction');
            polar(axe, plotx_dir(:)/180*pi, [data_ori(:); data_ori(1)]);
            set(axe ,'LineWidth',2);
        else
            plot(axe, self.para.condition_prop_list.direction, data_ori, '*b-','LineWidth', 1, 'MarkerSize', 1);
            ylim(lim)
        end
        
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