function current = removecurrents(data,samplerate,a0,a1)

% delete spikes from intracellular recording membrane potential trace
%     Input: ---- data,  was directly loaded from abf file, through abfload.m (pclamp 9) or import_abf.m(pclamp 10) 
%                       if abf file has multi-channels, then output will be  mp trace of channel 1, with spikes removed.
%            ---- samplerate,  default samplerate was 10k Hz
%            ---- a0, interpolation beginning point (ms),  default was 1 ms
%                     preceding spike time point, spike time point was determined by
%                     findspikes.m,  usually not the max point
%            ---- a1, interpolation end point (ms), default was 5 ms after spike time point
%     Output:---- membrane potential vector with spikes removed

% Yingjie Zhu.  2007/08/24  Modified
% Yingjie Zhu.  2007/11/15  in order to avoid error, extend mp to longer new_mp 


if nargin<2,samplerate = 10000;end; % default samplerate was 10k Hz
if nargin<3,a0=-1*samplerate/1000;end; % 1 ms befor spike time point 
if nargin<4,a1=5*samplerate/1000;end; % 5 ms after spike time point


mp = data(:,1); % in general, channel 1 was membrane potential trace.
mean_mp = mean(mp);
length_mp = length(mp);

%%% find spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spikes=findcurrents(mp,samplerate);
spike_time=find(spikes==1);
num_spikes=length(spike_time);

 
%%% Linear Interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to avoid error (index out of bounds), extend mp to new_mp by mean_mp
extend_mp = mean_mp*ones(a1,1);
new_mp = [mp;extend_mp];

for i=1:num_spikes    
ts=spike_time(i)+a0;
te=spike_time(i)+a1;
if ts <=0
    ts = 1;
end;
delta=(new_mp(te)-new_mp(ts))/(a1-a0);
q=a1-a0;
      for j=1:q
          new_mp(ts+j)=new_mp(ts)+j*delta;
      end;
end

current = new_mp(1:length_mp);

