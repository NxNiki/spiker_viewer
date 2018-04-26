% demo for data analysis based on class SV.

%  Get the orientation tuning curve from intracellular data.(according to spikes)
%   Input     ----- SV object;
%                   data file;(ch1: recording, ch2: stim tag, for selecting data of specific stimulus condition)
%   Output    ----- orientation, spatial, temporal frequency tuning curve.

clear
close all

input_abf = '/Users/Niki/Documents/MATLAB/Niki/spike_viewer/abf/20151118_cell1_1_0002_bino_fori_tv.abf';
input_mat = '/Users/Niki/Documents/MATLAB/Niki/spike_viewer/mat/001_Dur{3.000}Dir[12]Sf[3]Tf[2.000].mat';
output_dir = '/Users/Niki/Documents/MATLAB/Niki/spike_viewer/output/';
output_fig_ext = '.jpg';

dim.rep = 1;
dim.dir = 2;
dim.sf = 3;
dim.tf = 4;
dim.sp = 5;

[~,file_name] = fileparts(input_abf);

if ~exist(output_dir,'dir')
    mkdir(output_dir)
end

% Load ABF 
data = abfload(input_abf);

% get the stim parameter info:
load(input_mat, 'Para');
P = Para;clear Para
Para = para('onset', P.onset,'duration', P.value(1,1), 'direction', P.value(:,2), 'spatial_frequency',P.value(:,3), 'temporal_frequency', P.value(:,4));
SV = spike_viewer(Para, data(:,1), data(:,2), dim);
clear data

%--------------------------------------------------------------------------

SV.spike_thresh = 0.3;
SV.include_iti = true;

%--------------------------------------------------------------------------

mp_raw = SV.rec3d_repeat_para_sample;

condition_prop_len = SV.para.condition_prop_length;
num_tf = condition_prop_len(3);
num_sf = condition_prop_len(2);
num_dir = condition_prop_len(1);
ParaList = SV.para.condition_prop_list;
iti = SV.para.iti;
sr = SV.para.sample_rate;
duration = SV.para.duration;
plot_xtime = SV.plot_xtime;
unique_para_val = SV.para.unique_condition_prop_val;

%% fitting tuning curve:
SV.tuning_curve_bimodal_gaussian

clear
%% plot spike rate heat map:

SR = SV.spike_rate_hyper_array;
% y: direction, x: spatial frequency
for tf = 1:num_tf
    f = figure;
    % heatmap for each repetition
%     determine the position of plot:
%     [x,y]=subplotshape(num_rep);
%     set(gcf,'Color','w','position',[1, 1, y*300, x*225]);
%     for rep = 1:size(SR,dim.rep)
%         subplot(x,y,rep)
%         smooth_heatmap(squeeze(SR(rep,:,:,tf)), ParaList.sf, ParaList.dir);
%     end
    
    % mean heat map for all retetitions:
    %smooth_heatmap(squeeze(mean(SR(:,:,tf,:))), ParaList.sf, ParaList.dir);
    drawheatmap(squeeze(mean(SR(:,:,:,tf))), ParaList.spatial_frequency, ParaList.direction);
    colorbar
    suptitle(sprintf('temporal frequency: %0.2f Hz', ParaList.temporal_frequency(tf)))
    WritePlot(f, fullfile(output_dir, file_name, [sprintf('HeatMap_tf_%0.3f', tf), output_fig_ext]))
end

%% plot PSTH (peristimulus time histogram)

[ST] = SV.spike_trace_hyper_array;
psth_peak = nan(num_tf, num_sf, num_dir);
% pos = MakeSubplotPosition(num_dir*2, num_sf, 0.02, 0.03);
for i_tf = 1:num_tf
%     f = figure;
%     set(gcf,'Color','w','Position',[0 0 num_sf*400 num_dir*900]);
    for i_sf = 1:num_sf
        for i_dir = 1:num_dir
            
            f = figure;
            set(gcf,'Color','w','Position',[0 0 550 350]);
            d = find(squeeze(ST(:,i_dir,i_sf,i_tf,:))');
            
            %axe1 = subplot(2*num_dir, num_sf, ind, 'Position', pos(ind,:));
            axe1 = subplot(2, 1, 1);
            scatter(mod(d(:)-1,SV.sample_length)+1, ceil(d(:)/SV.sample_length),5,'filled')
            set(axe1,'Ylim',[0,SV.para.num_repeat+1])            
            set(axe1,'Ytick',1:SV.para.num_repeat)
            set(axe1,'Xlim', SV.plot_xlim)
            set(axe1,'Xtick',[])
            set(axe1,'XTicklabel',[])
            
            pos = get(axe1,'pos');
            trans_pos = [0, 0.1, 0, -0.1];
            set(axe1,'pos', pos+trans_pos)
            
            hold on
            plot([iti, iti]*sr, [0,SV.para.num_repeat+1],'r:');
            %hold on
            plot(([iti, iti] + duration)*sr, [0,SV.para.num_repeat+1],'r:');
            ylabel(axe1, 'Trials')
            
            %axe2 = subplot(2*num_dir, num_sf, ind + num_sf, 'Position', pos(ind+num_sf,:));
            axe2 = subplot(2, 1, 2);
            [~, r] = psth(axe2, d(:), 100, SV.para.sample_rate, SV.para.num_repeat, SV.sample_length, true);
            psth_peak(i_tf, i_sf, i_dir) = max(r);
            y_lim = [-5, 100];
            set(axe2, 'Ylim', y_lim);
            hold on
            plot([iti, iti]*sr/10, y_lim,'r:');
            %hold on
            plot(([iti, iti] + duration)*sr/10, y_lim,'r:');
            
            set(axe2,'Xtick',SV.plot_xtick/10)
            set(axe2,'XTicklabel',SV.plot_xticklabel)
            
            pos = get(axe2,'pos');
            trans_pos = [0, 0, 0, .2];
            set(axe2,'pos', pos+trans_pos)

%             if i_dir == num_dir
                xlabel(axe2, sprintf('spatial frequency: %0.3f', ParaList.spatial_frequency(i_sf)))
%             end
%             if i_sf == 1
                ylabel(axe2, sprintf('direction: %3d', ParaList.direction(i_dir)))
%             end

            ind = (i_dir-1)*2*num_sf + i_sf;
            WritePlot(f, fullfile(output_dir, file_name, [sprintf('PSTH_tf_%0.3f_%02d', i_tf, ind), output_fig_ext]))
        end
    end
%     suptitle(sprintf('temporal frequency: %1.3f', ParaList.tf(i_tf)));
%     WritePlot(f, fullfile(output_dir, file_name, [sprintf('PSTH_tf_%0.3f', i_tf), output_fig_ext]))
end              


%% plot raw trace, f1
% determine the position of plot:
[x,y]=subplotshape(SV.para.num_unique_stm);
figure;
set(gcf,'Color','w','Position',[0 0 x*300 y*225]);
hst = suptitle(file_name);
set(hst,'Position',[0.549267 -0.0127427 9.16025],'Interpreter','none');

color = {'r','c','m','g'}; %red, cyan, magenta, green
for i = 1:SV.para.num_unique_stm
    subplot(x,y,i);
    
    mp_trace = squeeze(mp_raw(:,i,:));
    % shift the trace in order to plot all repeats in the same subplot
    mp_trace_shift = bsxfun(@plus, mp_trace, (0:SV.para.num_repeat-1)'*2);
    
    for j = 1:SV.para.num_repeat
        
        hp = plot(mp_trace_shift(j,:));
        hold on;
        if j>1
            set(hp,'Color',color{j-1}); 
        end
    end
    set(gca,'Box','off','YLim',[min(mp_trace_shift(:))-1 max(mp_trace_shift(:))+1]);
    title(sprintf('ori:%d, sf: %0.3f, tf %0.3f', unique_para_val(i,:)));
end

if ispc
    saveppt(fullfile(output_dir, [file_name,'.ppt']));
end

%% plot average trace, f2
mp_trace_ave = squeeze(mean(mp_raw));
% determine the position of plot:
[x,y]=subplotshape(SV.para.num_unique_stm);
figure;
set(gcf,'Color','w','position',[1, 1,x*300, y*225]);

for j=1:SV.para.num_unique_stm
    subplot(x,y,j);
    plot(plot_xtime, mp_trace_ave(j,:));
    hold on;
    plot([round(iti), round(iti)]/sr, [min(mp_trace_ave(:)) max(mp_trace_ave(:))],'r:');
    hold on;
    plot([round(iti), round(iti)]/sr + duration, [min(mp_trace_ave(:)) max(mp_trace_ave(:))],'r:');
    set(gca,'Box','off');
    xlim([min(plot_xtime), max(plot_xtime)]);
    xlabel('sec');ylabel('mv');
    title(strcat('direction = ',num2str(unique_para_val(j))));
end;
suptitle('average mp trace');

if ispc
    saveppt(fullfile(output_dir, [file_name,'.ppt']));
end

%% plot tuning curve, f3
% polar plot with spikes

f = SV.tuning_curve('spike', 'polar');
WritePlot(f, fullfile(output_dir, file_name, ['turning_curve_sr_polarplot', output_fig_ext]))

%% polar plot with mp

f = SV.tuning_curve('m.p.','polar');
WritePlot(f, fullfile(output_dir, file_name, ['turning_curve_mp_polarplot', output_fig_ext]))

%% calibrate the oritation selectivity index with spike

f = SV.tuning_curve('spike');
WritePlot(f, fullfile(output_dir, file_name, ['turning_curve_sr', output_fig_ext]))

%% calibrate the oritation selectivity index with membrane potential

f = SV.tuning_curve('spike');
h_tmp = suptitle(file_name);
set(h_tmp,'Interpreter','none');
WritePlot(f, fullfile(output_dir, file_name, ['turning_curve_mp', output_fig_ext]))

%% get parameter and corresponding OSI:
[osi, para] = SV.osi_report('spike');

%% time frequency analysis:
[f,Y,isSimple,ParaVal] = SV.time_frequency;
num_plot = size(Y,1);
[plotx,ploty] = subplotshape(num_plot);
y_lim = [min(Y(:)), max(Y(:))];
figure
set(gcf,'Color','w','position',[1, 1,plotx*300, ploty*225]);
for i=1:num_plot
    subplot(plotx,ploty,i)
    plot(f,Y(i,:))
    ylim(y_lim)
    title(sprintf('dir: %d,sf:%1.2f, tf:%1.0f',ParaVal(i,1),ParaVal(i,2),ParaVal(i,3)))
end

%% draw spike shape:
[all_spike_shape, average_spike_shape, spike_dur, plotx] = SV.spike_shape(.2,0);
[spike_range, peak_valley_ratio, peak_valley_duration, recoverSlope, valleyFWHM] = SV.spike_shape_analysis;
figure
plot(plotx, all_spike_shape,'Color',[.75,.75,.75])
hold on
plot(plotx, average_spike_shape,'Color',[0,0,0],'LineWidth',1)
set(gca,'YAxisLocation','Origin')
set(gca,'XAxisLocation','Origin')
set(gca,'fontsize',14);

set(gca,'Box','Off')
plot([10; 10], [1; 1.5], '-k',  [10; 15], [1; 1], '-k', 'LineWidth', 2)
text(10,1.25, '.5 mV ', 'HorizontalAlignment','right','FontSize',15)
text(12.5,0.95, '.5 ms', 'VerticalAlignment','top','HorizontalAlignment','center','FontSize',15)
% set(gca, 'Visible', 'off')
set(gcf,'Color','w','position',[1, 1, 550, 550]);
hold off
%% save

% if ispc
%     xlswrite(strcat(output_dir, filesep, file_name,'OSI_spike'),{file_name,OSI_spike},1);
% end


