function NaCurrents = findcurrents(currents,samplerate,Thresh,MaxCurDur)

% FINDCURRENTS finds the sodium current times in voltage-clamped traces
%
%	currenttimes = findcurrents(currents,samplerate,Thresh,MaxCurDur)
%
%	finds the sodium currents by high-pass filtering the currents and 
%	finding the crossings of the result with a threshold Thresh (in pA).
%	The crossings can be at most MaxCurDur ms apart to qualify. 
%	It returns currents, a sparse matrix which is 1 where it found a inward sodium current.
%
%	currens is a matrix nsamples by nstimuli, in pA
%	samplerate is in Hz
%	Thresh is in pA 		(default: -200 pA)
%	MaxSpkDur is in ms 	(default: 4 ms)
%
% 1996 David Ferster and Matteo Carandini
% 1997 Matteo Carandini
% part of the Matteobox toolbox
% 
% modified by Xiaojing Aug.18,2008

if nargin < 4, MaxCurDur = 4; end
if nargin < 3, Thresh = -200; end

[nsamples, nstimuli] = size(currents);

% --------------- define the highpass filter 
cutoff = 1000;	% Hz. Not really. With 1000 it cuts off around 100 Hz.

% 1 - two RC circuits in series
deltat = 1/samplerate;
tau = 1/cutoff;
b = [ 1-(deltat/tau)^2 2*(deltat-tau)/tau ((deltat-tau)/tau)^2];
a = [ 1 2*(deltat-tau)/tau ((deltat-tau)/tau)^2];

% to look at the filter: freqz(b,a);

% -------------------- high-pass filter the potentials:
hipasdpots = zeros(nsamples,nstimuli);
% hipasdpots(:) = filter( b,a, potentials(:) ); 
% do it one by one. Filter is memory intensive.
for istim = 1:nstimuli
	hipasdpots(:,istim) = filter( b,a, currents(:,istim) ); 
end

% ------------------- if you want to see the residuals
% figure; plot(potentials(:,1)-hipasdpots(:,1))
% figure; plot(hipasdpots(:,1))

%-------------------- find the spikes
NaCurrents = sparse([],[],[],nsamples,nstimuli,1000); %% at most 1000 spikes/stim
for istim = 1:nstimuli
	above = find( hipasdpots(2:nsamples-1,istim)<Thresh );
	if any(above)
		diffabove = diff(above);
		spbegs = above([ 1; find(diffabove>1)+1 ]);
		spends = above([ find(diffabove>1); length(above) ])+1;
		sptimes = spbegs((spends-spbegs)<(MaxCurDur*samplerate/1000));
		if any(sptimes)
			NaCurrents(sptimes,istim) = ones(size(sptimes));
		end
	end
end

%-------------------- look at the results:
% clf;
% istim = 7;
% plot (hipasdpots(:,istim));
% hold on;
% plot ( find(spikes(:,istim)),hipasdpots(spikes(:,istim),istim), 'go');



