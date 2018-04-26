classdef para
    %PARA Para store parameters of stimuli that are necessary for data
    %analysis, which is implemented by spikeviewer. 
    %usage:
    %   P = para('onset', Onset_val, 'spatial_freq', SF_val, ...)
    %   by Niki 2015/7/18.
    
    properties
        sample_rate = 10^4; % in Hz
        pre_iti = 0.25;
        pos_iti = 0.25;
    end
    
    properties (SetAccess = private)
        onset
        duration
        condition_prop
        
        iti % para set iti as the minimal value in pre- and post-iti
        num_stm
        num_unique_stm
        num_repeat
        num_condition_prop
        condition_prop_name
        condition_prop_val
        condition_prop_length
        condition_prop_list
        sort_trials_ind
        unique_condition_prop_val
    end
    
    properties(Dependent = true, SetAccess = private)
        unique_stm_ind % num_unique_stm * num_repeat matrix, contains index of stm coresponding to unique_stm_value
        repeat_ind
        trial_length
        iti_length % samples of iti, same for pre and post iti (min value).
    end
    
    methods
        function P = para(varargin)
            % assgin input to properties:
            num_input_prop = nargin/2;
            if mod(num_input_prop,1)~=0
                error('PARA:input parameters should be pairs of name and value')
            else
                for i=1:num_input_prop
                    % make sure every property is a column vector or scalar:
                    prop_name = varargin{i*2-1};
                    prop_value = varargin{i*2};
                    try
                        P.(prop_name) = prop_value(:);
                    catch
                        disp(['PARA:', prop_name, 'is not a property of class para. Adding it to condition properties'])
                        P.condition_prop.(prop_name) = prop_value(:);
                    end
                end
            end
            
            P.condition_prop_name = fieldnames(P.condition_prop);
            % check if some parameters have equal length of onset:
            P.num_stm = length(P.onset);
            if P.num_stm<=0
                error('PARA:property onset not assigned')
            end
            
            condition_prop_val = cell2mat(struct2cell(P.condition_prop)');
            
            % check if length of condition property is same
            P.num_condition_prop = length(P.condition_prop_name);
            
            [P.condition_prop_val, P.sort_trials_ind] = sortrows(condition_prop_val, P.num_condition_prop:-1:1);
            
            [P.unique_condition_prop_val] = sortrows(unique(P.condition_prop_val, 'rows'), P.num_condition_prop:-1:1);
            P.num_unique_stm = size(P.unique_condition_prop_val,1);
            P.num_repeat = P.num_stm/P.num_unique_stm;
            
            if mod(P.num_repeat,1)~=0
                error('PARA: num_repeat is not a round value, indicating each condition do not have same number of trials')
            end
            
            P.condition_prop_list = structfun(@unique, P.condition_prop, 'UniformOutput', false); % return cell array
            P.condition_prop_length = structfun(@length, P.condition_prop_list); % return vector
            
            if ~isempty(P.pre_iti)||~isempty(P.pos_iti)
                P.iti = min([P.pre_iti(:); P.pos_iti(:)]);
            else
                P.iti = 0.1;
            end
            
        end
        
        % set methods:-----------------------------------------------------
        
        function self = set.onset(self, ons)
            if any(diff(ons)<=0)
                error('PARA:onset should be vector of increasing floating numbers')
            end
            
            self.onset = ons;
        end
        
        % get methods:-----------------------------------------------------
        
        function out = get.trial_length(self)
            out = self.duration * self.sample_rate;
        end
        
        function out = get.iti_length(self)
            out = self.iti * self.sample_rate;
        end
        
        function out = get.repeat_ind(self)
            out = find_rep_ind(self.condition_prop_val);
        end
        
        % -----------------------------------------------------------------
        function out = find_condition_prop(self, para_line, option)
            if size(para_line,2)~=self.num_condition_prop
                error('PARA: input vector should have same number of columns as number of condtion properties')
            end
            
            if nargin <3
                option = 'all';
            end
            
            % regard nan columns as selecting all values:
            col = ~isnan(para_line)&~isempty(para_line);
            
            switch option
                case 'all'
                    out = ismember(self.condition_prop_val(:,col), para_line(col),'rows');
                case 'unique'
                    out = ismember(self.unique_condition_prop_val(:,col), para_line(col), 'rows');
            end
        end
        
    end
    
end

