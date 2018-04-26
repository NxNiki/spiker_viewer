function [num_figures, num_rows, num_cols] = plot_position(self, condition_prop_row, condition_prop_col)

if nargin < 2
    % plot all axes in one figure
    num_figures = 1;
    [num_rows, num_cols] = subplotshape(self.num_unique_stm);
    return
end

if nargin == 3
    all_condition_props = self.para.condition_prop_name;
    prop_across_figure = setdiff(all_condition_props, [condition_prop_row, condition_prop_col]);
    num_figure = self.para.(prop_across_figure)