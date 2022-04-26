% Testing a grid cell CAN network with IZ neurons
% Nate Sutton 2022
simdur = 1050;%1010;%490;%1300;%100e3; % total simulation time, ms
spiking_bin = 10;

ncells = Ne; % total number of cells per layer
Ne=ncells; Ni=ncells;
a=[0.02*ones(Ne,1)];
b=[0.2*ones(Ne,1)];
c=[-65*ones(Ne,1)];
d=[8*ones(Ne,1)];
ai=[0.02*ones(Ni,1)];
bi=[0.25*ones(Ni,1)];
ci=[-65*ones(Ni,1)];
di=[2*ones(Ni,1)];
p = [a, b, c, d];
p2 = [ai, bi, ci, di];
tau = 30; %% Cell parameters % grid cell synapse time constant, ms
t = 0; % simulation time variable, ms
skip_t = 10; % initial time to skip because pregenerated initial firing is loaded in this time
v=-65*ones(Ne,1); % Initial values of v
u=b.*v;
vi=-65*ones(Ni,1); % inhib neurons
ui=b.*vi;
load('Ii_initial2.mat'); % initial gc firing
Ii = Ii_initial2;
load('init_firings3.mat'); % initial gc firing
firings = init_firings3;
load('../../move_test3/data/B_saved.mat'); % velocity input matrix
Ie=60*(B.^1.8)'; % excitatory input
load('../../move_test3/data/mex_hat3.mat'); % load weight matrix
mex_hat = mex_hat3*2.8;
mex_hat = mex_hat-0.0015;
mex_hat = mex_hat.*(mex_hat>0); % no negative values
% video parameters
ccol = load('../../move_test3/src/neuron_space_colormap.mat');
savevideo = true;
h = figure('color','w','name','');
%numberOfFrames = simdur-skip_t;
numberOfFrames = (simdur-skip_t)/spiking_bin;
allTheFrames = cell(numberOfFrames,1);
vidHeight = 337;%342;
vidWidth = 442;%434;
allTheFrames(:) = {zeros(vidHeight, vidWidth, 3, 'uint8')};
allTheColorMaps = cell(numberOfFrames,1);
allTheColorMaps(:) = {zeros(256, 3)};
myMovie = struct('cdata', allTheFrames, 'colormap', allTheColorMaps);
set(gcf, 'nextplot', 'replacechildren'); 
set(gcf, 'renderer', 'zbuffer');
caxis manual;          % allow subsequent plots to use the same color limits

for t=skip_t:simdur % simulation of 1000 ms
	fired=find(v>=30); % indices of spikes
	firings=[firings; t+0*fired,fired];
	Ii = inhib_curr(Ii, t, mex_hat, firings, tau);
	%Ie = Ie .* (1 + (rand(ncells,1)*.02)); % add some random noise
	%Ie = Ie .* (1 - (rand(ncells,1)*.02)); % subtract some random noise
	[v, u] = iznrn(v, u, p, fired, Ie, Ii);
	if savevideo & mod(t,spiking_bin) == 0
		myMovie = heatmap(ncells, firings, t, skip_t, h, myMovie, ccol, spiking_bin);
    end
end

close(h);
if savevideo
	videofile = VideoWriter('heatmap.avi'); % Create a VideoWriter object to write the video out to a new, different file.
	open(videofile)
    myMovie(1) = [];
	writeVideo(videofile,myMovie) % Write the movie object to a new video file.
	close(videofile)
end

function [v, u] = iznrn(v, u, p, fired, Ie, Ii)
	a=p(:,1);b=p(:,2);c=p(:,3);d=p(:,4);
	v(fired)=c(fired);
	u(fired)=u(fired)+d(fired);
	v=v+(0.04*v.^2+5*v+140-u+Ie-Ii); % step 1.0 ms
	u=u+a.*(b.*v-u);
end

function spikes = fbin(ni, startt, endt, firings) 
	% report firing in time bin
	% ni = neuron index, startt = start time, endt = end time
	st = []; % spike times
	all_firing = (find(firings(:,2)==ni));
	for si=1:size(all_firing)
		st = [st; firings(all_firing(si),1)];
	end
	st = find(st>startt & st<endt);
	spikes = size(st);
end

function spike_found = find_spike(ni, t, firings) 
	% report spike times in a bin of time
	spike_found = false;
	all_spike_times = (find(firings(:,2)==ni));
	for si=1:size(all_spike_times)
		spike_time = firings(all_spike_times(si),1);
        if spike_time == t
            spike_found = true;
        end
    end
end

function Ii = inhib_curr(Ii, t, mex_hat, firings, tau)
	% generate inhibitory currents
	gc_firing = zeros(size(mex_hat,1)); 
	for i=1:size(Ii)
        spike_found = find_spike(i,t,firings);
		if spike_found == true
	        gc_firing(:,i) = gc_firing(:,i)+1;
	    end    
    end
    in_current = (mex_hat*gc_firing')';

	% calculate tau factor
	o = ones(size(mex_hat(:,1)));
	in_summed = in_current'*o; in_summed2 = in_summed;
	in_summed2 = in_summed2.*(in_summed2>0); % no negative values
	Ii = Ii + (in_summed2 - Ii)/tau;

	% todo: replacing gc_firing*cen_sur=in_current_to_gc
	% gc_firing*syn_wt=current_to_in; iz_nrn(current_to_in)=in_firing;
	% in_firing*cen_sur=current_to_gc
end

function myMovie = heatmap(ncells, firings, t, skip_t, h, myMovie, ccol, spiking_bin)
	binned_firing = [];
	for i=1:sqrt(ncells)
		temp = [];
		for j=1:sqrt(ncells)
			nrn_i = ((i-1)*sqrt(ncells))+j;
			spk_t = fbin(nrn_i, t-10, t, firings);
			temp = [temp; spk_t(1,1)];
		end
		binned_firing = [binned_firing; temp'];
	end
	%figure(h);
	cla reset;
  	hAxes = gca;
	imagesc(hAxes, binned_firing);
	colormap(ccol.CustomColormap2);
	xlabel('neuron position on x axis') 
	ylabel('neuron position on y axis')
	shading interp;
	%axis square
    axis('tight')
	title({sprintf('t = %.1f ms',t),'Population activity'})
	set(gca,'ydir','normal')
    caxis([0 8.0])
	cb = colorbar;
	set(cb, 'ylim', [0 5.5]); % set colorbar range
	drawnow
	thisFrame = getframe(gcf);
    %frames_number = (t-(skip_t-1))/spiking_bin + 1
    frames_number = t/spiking_bin;
  	myMovie(frames_number) = thisFrame;
end