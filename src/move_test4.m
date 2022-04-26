% Testing a grid cell CAN network with IZ neurons
% Nate Sutton 2022
clear all;
clc;
simdur = 500;%1010;%490;%1300;%100e3; % total simulation time, ms
spiking_bin = 10;

ncells = 900; % total number of cells per layer
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
tau = 35; %% Cell parameters % grid cell synapse time constant, ms
gcintau = 35;
ingctau = 35;
t = 0; % simulation time variable, ms
skip_t = 10; % initial time to skip because pregenerated initial firing is loaded in this time
v=-65*ones(Ne,1); % Initial values of v
u=b.*v;
vi=-65*ones(Ni,1); % inhib neurons
ui=bi.*vi;
load('Ii_initial2.mat'); % initial gc firing
Ii = Ii_initial2;
Ie2 = zeros(ncells,1);
Ii3 = zeros(ncells,1);
load('init_firings4.mat'); % initial gc firing
firings = init_firings4;
gc_firings = [10,5];
in_firings = [10,5];
load('../../move_test3/data/B_saved.mat'); % velocity input matrix
Ie=68*(B.^10)'; % excitatory input
load('../../move_test3/data/mex_hat3.mat'); % load weight matrix
mex_hat = mex_hat3*3;
mex_hat = mex_hat-0.0022;
mex_hat = mex_hat.*(mex_hat>0); % no negative values

if false
    max_ind = sqrt(size(mex_hat(:,1)));
    for x = 1:max_ind
        for y = 1:max_ind
            ind = ((y-1) * max_ind) + x;
            if (get_pd(x,y) ~= 'r')                
                %mex_hat(:,ind) = .1;
                %mex_hat(:,ind) = 0;
                Ie(ind) = 1.1;%1.0;
                %B(ind) = .5;
            else
                Ie(ind) = 1.0;%.5;
            end
            %B(ind) = 1;
        end
    end
end
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
	gc_fired=find(vi>=30); % indices of spikes
	gc_firings=[gc_firings; t+0*gc_fired,gc_fired];
    [Ie2, vi, ui] = gc_in_signal(Ie2, t, gc_fired, gc_firings, vi, ui, p2, gcintau, ncells);
	in_fired=find(vi>=30); % indices of spikes
    in_firings=[in_firings; t+0*in_fired,in_fired];
    Ii3 = in_gc_signal(t, mex_hat, in_firings, ncells, Ii3, gcintau);

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

function [Ie2, vi, ui] = gc_in_signal(Ie2, t, gc_fired, gc_firings, vi, ui, p2, gcintau, ncells)
    % generate inhibitory currents
    o = ones(ncells,1);    
	gc_firing = zeros(ncells); 
	for i=1:ncells
        spike_found = find_spike(i,t,gc_firings);
		if spike_found == true
	        gc_firing(:,i) = gc_firing(:,i)+1;
	    end    
    end
	% todo: replacing gc_firing*cen_sur=in_current_to_gc
	% gc_firing*syn_wt=current_to_in; iz_nrn(current_to_in)=in_firing;
	% in_firing*cen_sur=current_to_gc
    gc_to_in_wt = 1; % gc to in synapse weight
	gcin_current = (gc_to_in_wt*gc_firing')';
	gcin_summed = gcin_current'*o;    
    Ie2 = Ie2 + (gcin_summed - Ie2)/gcintau;
	Ii2 = zeros(1,900)';
	[vi, ui] = iznrn(vi, ui, p2, gc_fired, Ie2, Ii2);
end

function Ii3 = in_gc_signal(t, mex_hat, in_firings, ncells, Ii3, gcintau);
	% generate in to gc signalings
	o = ones(ncells,1);
	in_firing = zeros(ncells); 
	for i=1:ncells
        spike_found = find_spike(i,t,in_firings);
		if spike_found == true
	        in_firing(:,i) = in_firing(:,i)+1;
	    end    
    end
    in_current = (mex_hat*in_firing')';
    in_to_gc_wt = 1; % in to gc synapse weight
	ingc_current = ((mex_hat*in_to_gc_wt)*in_firing')';
	ingc_summed = ingc_current'*o;  
    Ii3 = Ii3 + (ingc_summed - Ii3)/gcintau;
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
    caxis([0 4.0])
	cb = colorbar;
	set(cb, 'ylim', [0 5.5]); % set colorbar range
	drawnow
	thisFrame = getframe(gcf);
    %frames_number = (t-(skip_t-1))/spiking_bin + 1
    frames_number = t/spiking_bin;
  	myMovie(frames_number) = thisFrame;
end

function pd = get_pd(x, y)
    % find neuron preferred direction
	if (mod(y,2) == 0)
		if (mod(x,2) == 0)
			pd = 'd';
		else 
			pd = 'r';
        end
    else
		if (mod(x,2) == 0)
			pd = 'l';
        else
			pd = 'u';	
        end
    end
end