classdef spike_viewer
    
    %SPIKE_VIEWER extract spikes and membrane potenatials corresponding to each
    %unique combination of stimulus parameters.    
    % INPUT: para ojbect containing information of stimuli parameters.
    % use this class as: out = stim_info(para);
    % see also para.m    
    % by Niki 2015/6/29.
    
    % make include_iti as a properties to avoid repetitious use in methods.
    % by Niki 2015/11/7.    
    % add propertiy spontaneous_spike_rate; rec2d_para_sample return
    % spontaneous_rec.
    % by Niki 2016/03/18.
    
    properties
        
        stm_thr_rat = 2/3; % used to determine onset sample.
        refractory_period = 4; % in milliseconds
        include_iti = false; % if true, rec2d and rec3d sample will include pre and post iti.
        spike_thresh = 1; % parameters for function findspikes
        osi_type = 2; % [1,2] two way of computing osi.
        dim_info
        num_bins = 20; % number of bins in psth plot
        spike_rate_minus_spontaneuous = true; %
        spike_peak_alignment = true;
        
    end
    
    properties (SetAccess = private)
        
        para
        recording % N*1 vector of neural response. For recording system in Xiaohui Zhang's lab, this should be the 1st column of recording data.
        stml % same size as recording, used to determine the onset tag. For recording system in Xiaohui Zhang's lab, this should be the 2nd column of recording data.
        dim_val
        
    end
    
    properties (Dependent = true, SetAccess = private)
        
        sample_length % number of samples for single stimulus. (including pre and post iti length)
        sample_index % if includeiti, the index will start from negative integer, else, start from 0.
        sample_duration
        start_sample_index
        spontaneous_recording
        spontaneous_spike_trace
        spontaneous_spike_rate% spike rate before the display of stimuli, use same duration as sample_duration
        
    end
    
    methods
        
        % constructing methods:--------------------------------------------
        
        function SV = spike_viewer(p, rec, stm, dif)
            
            SV.para = p;
            
            SV.recording = rec;
            
            SV.stml = stm;
            
            SV.dim_val = cat(2, {1:p.num_repeat}, p.condition_prop_list, {SV.sample_index});
            
            SV.dim_info = dif;
            
        end
        
        % get methods:-----------------------------------------------------
        
        function out = get.sample_length(self)
            
            out = self.para.trial_length;
            if self.include_iti
                out = out + self.para.iti_length * 2;
            end
            
        end
        
        function out = get.sample_index(self)
            
            out = 0:self.sample_length-1;
            if self.include_iti
                out = out - self.para.iti_length;
            end
            
        end
        
        function out = get.sample_duration(self)
            
            out = self.sample_length/self.para.sample_rate;
            
        end
        
        function out = get.start_sample_index(self)
            
            thresh = self.stm_thr_rat * max(self.stml);
            
            out = find(self.stml>thresh,1,'first');
            
        end
        
        function out = get.spontaneous_recording(self)
            
            start_sample = self.start_sample_index;
            
            spl = self.sample_length;
            
            out = self.recording(start_sample - spl: start_sample-1);
            
        end
        
        function out = get.spontaneous_spike_trace(self)
            
            spon_rec = self.spontaneous_recording;
            
            out = findspikes(spon_rec(:), self.para.sample_rate, self.spike_thresh, self.refractory_period);
            
        end
        
        function out = get.spontaneous_spike_rate(self)
            
            out = sum(self.spontaneous_spike_trace)/self.sample_duration;
            
        end
        
        % others:----------------------------------------------------------
        
        function self = add_parameter(self, varargin)
            
            for i = 2:2: nargin
                
                self = self.add_para(varargin{i},varargin{i+1});
                
            end
            
        end
        
        function out = extract_recording(self, para_line)
            
            ind = self.para.find_row_ind(para_line);
            
            rec2d = self.rec2d_para_sample;
            
            out = rec2d(ind,:);
            
        end
        
        function rec2d = rec2d_para_sample(self)
            
            % convert 1-D recordings to 2-D matrix. Each row corresponds to
            % one single stimulus.
            % index (relative to start_sample) of sample corresponding to the onset of each single stimulus
            
            onset_sample = round(self.para.onset * self.para.sample_rate) + self.start_sample_index;
            
            if self.include_iti
                
                onset_sample = onset_sample - self.para.iti_length;
                
            end
            
            spl = self.sample_length;
            
            tag_ind = bsxfun(@plus,repmat(onset_sample,[1,spl]), 0:spl-1);
            
            rec2d = reshape(self.recording(tag_ind(:)), self.para.num_stm, spl);
            
            % make recording the same order as para: this step is crucial
            % for otherwise parameters will not correspond with recording
            % correctly
            
            rec2d = rec2d(self.para.sort_trials_ind,:);
            
        end
        
        function rec3d = rec3d_repeat_para_sample(self)
            
            % convert 2-d recording to 3-d array, with dimension 1:
            % repetition, dimension 2: trial index correponding to
            % para.unique_condition_prop_val, dimension 3: sample length
            
            rec2d = self.rec2d_para_sample;
            
            sp_length = size(rec2d,2);
            
            rec3d = reshape(rec2d, self.para.num_repeat, self.para.num_unique_stm, sp_length);
            
        end
        
        function out = membrane_potential(self)
            
            rec3d = self.rec3d_repeat_para_sample;
            
            sp_length = self.sample_length;
            
            mp = removespikes(rec3d(:));
            
            mp_trace = reshape(mp, self.para.num_repeat, self.para.num_unique_stm, sp_length);
            
            out = mean(mp_trace,3);
            
        end
        
        function out = spike_trace(self)
            
            % out is a num_repeat * num_unique_stm *sample_length matrix.
            
            rec3d = permute(self.rec3d_repeat_para_sample,[3,1,2]);
            
            sp_length = size(rec3d,1);
            
            out = findspikes(rec3d(:), self.para.sample_rate, self.spike_thresh, self.refractory_period);
            
            out = full(out); % convert to full array as ND sparse arrays are not supported.
            
            out = reshape(out, sp_length, self.para.num_repeat, self.para.num_unique_stm);
            
            out = permute(out, [2,3,1]);
            
        end
        
        function out = spike_rate(self)
            
            % out is a num_repeat * num_unique_stm matrix.
            
            self.include_iti = false;
            
            spike_trace = self.spike_trace;
            
            if self.spike_rate_minus_spontaneuous
                
                out = sum(spike_trace,3)/self.para.duration - self.spontaneous_spike_rate;
                
            else
                    
                out = sum(spike_trace,3)/self.para.duration;  
                    
            end
            
        end
        
        function out = spike_trace_hyper_array(self)
            
            % put condition properties (e.g. spatial frequency, temporal
            % freqency and direction) in separate dimensions, so that data
            % for each combination of condition can be easily selected.
            
            st = self.spike_trace;
            
            dim = self.para.condition_prop_length;
            
            out = reshape(st, [self.para.num_repeat, dim(:)', self.sample_length]);
            
        end
        
        function out = spike_rate_hyper_array(self)
            
            sr = self.spike_rate;
            
            dim = self.para.condition_prop_length;
            
            out = reshape(sr, [self.para.num_repeat, dim(:)']);
            
        end
        
        function out = membrane_potential_hyper_array(self)
            
            mp = self.membrane_potential;
            
            dim = self.para.condition_prop_length;
            
            out = reshape(mp, [self.para.num_repeat, dim(:)']);
            
        end

    end
    
end




