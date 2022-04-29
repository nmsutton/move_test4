% Testing a grid cell CAN network with IZ neurons
% Nate Sutton 2022
clear all;
clc;
simdur = 2010;%530;%3010;%210;%130;%1010;%490;%1300;%100e3; % total simulation time, ms
spiking_bin = 40;%40;

ncells = 900; % total number of cells per layer
Ne=ncells; Ni=ncells;
a=[0.1*ones(Ne,1)];%a=[0.02*ones(Ne,1)];
b=[0.2*ones(Ne,1)];
c=[-65*ones(Ne,1)];
d=[8*ones(Ne,1)];
ai=[0.1*ones(Ni,1)];%ai=[0.02*ones(Ni,1)];
bi=[0.25*ones(Ni,1)];
ci=[-65*ones(Ni,1)];
di=[2*ones(Ni,1)];
p = [a, b, c, d];
p2 = [ai, bi, ci, di];
tau = 35; %% Cell parameters % grid cell synapse time constant, ms
gcintau = 30;%35;
ingctau = 30;%35;
t = 0; % simulation time variable, ms
skip_t = 10; % initial time to skip because pregenerated initial firing is loaded in this time
v=-65*ones(Ne,1); % Initial values of v
u=b.*v;
vi=-65*ones(Ni,1); % inhib neurons
ui=bi.*vi;
load('Ii_initial2.mat'); % initial gc firing
Ii = Ii_initial2;
load('gc_ie_initial.mat'); % initial gc firing
gc_ie = gc_ie_initial;
%gc_ie = Ii_initial2;%zeros(ncells,1);
load('in_ii_initial.mat'); % initial gc firing
in_ii = in_ii_initial;
%in_ii = zeros(ncells,1);
load('init_firings4.mat'); % initial gc firing
firings = init_firings4;
gc_firings = init_firings4;%[10,5];
load('in_firings_init.mat');
in_firings = in_firings_init;
load('../../move_test3/data/B_saved.mat'); % velocity input matrix
%ext_ie=68*(B.^20)'; % excitatory input
ext_ie=60*(B.^22)'; % excitatory input
pd_match=68;%67;%132;%58;%%115;%89;%88;%74;%71.5;%83;%83;%60;%51;%42;%34.4;%43;
pd_nonmatch=60;%90;%30;%80;%60;
load('../../move_test3/data/mex_hat3.mat'); % load weight matrix
mex_hat = mex_hat3*3;
mex_hat = mex_hat-0.0022;
mex_hat = mex_hat.*(mex_hat>0); % no negative values
gc_to_in_wt = 25;%180;%25;%36;%47;%100;%180;%180;%30;%39;%180;%0.4;%0.2;%0.121;%;//0.12;%0.15; % gc to in synapse weight
in_to_gc_wt = 70;%410;%1200;%410;%410;%.45;%.45;%.39;%.15;%.15;%.3;%.15; % in to gc synapse weight

% tm model synapse parameters
global cap_ue tau_ue tau_xe tau_de gei u_ei x_ei ...
	   cap_ui tau_ui tau_xi tau_di gie u_ie x_ie;
cap_ue = .7;%.8;%9;%0.2; % U, utilization
tau_ue = 15;%40.0; % U signal decay time constant
tau_xe = 7.5;%15;%30;%100.0; % x signal decay time constant
tau_de = 30.0; % x signal decay time constant
gei = 1.0;
cap_ui = .8;%1;%.8;%9;%0.2; % U, utilization
tau_ui = 50;%30;%40.0; % U signal decay time constant
tau_xi = 30;%30;%60;%30;%15;%30;%100.0; % x signal decay time constant
tau_di = 40.0; % x signal decay time constant
gie = 1.0;
u_ei = zeros(ncells,1); % u before spike update
x_ei = ones(ncells,1); % x before spike update
u_ie = zeros(ncells,1); % u before spike update
x_ie = ones(ncells,1); % x before spike update
%i = 1;
%A = 2.0;

if true
    max_ind = sqrt(size(mex_hat(:,1),1));
    for x = 1:max_ind
        for y = 1:max_ind
            ind = ((y-1) * max_ind) + x;
            if (get_pd(x,y) == 'r')       
                ext_ie(ind) = pd_match;               
            else
                ext_ie(ind) = pd_nonmatch;
            end
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
gc_voltage=[]; in_voltage=[];
nrn_monit=702;%435;%772;%702;
in_fired=[];

for t=skip_t:simdur % simulation of 1000 ms
	gc_fired=find(v>=30); % indices of spikes
	gc_firings=[gc_firings; t+0*gc_fired,gc_fired];
    [gc_ie, vi, ui] = gc_in_signal(gc_ie, t, gc_firings, in_fired, vi, ui, p2, gcintau, ncells, nrn_monit, gc_to_in_wt);
	in_fired=find(vi>=30); % indices of spikes
    in_firings=[in_firings; t+0*in_fired,in_fired];
    [in_ii, in_firings] = in_gc_signal(t, mex_hat, in_firings, ncells, in_ii, gcintau, nrn_monit, in_to_gc_wt);
    [v, u] = iznrn(v, u, p, gc_fired, ext_ie, in_ii);
    in_voltage(end+1)=vi(nrn_monit);

	%Ii = inhib_curr(Ii, t, mex_hat, gc_firings, tau, nrn_monit);
	%[v, u] = iznrn(v, u, p, gc_fired, ext_ie, Ii);
	gc_voltage(end+1)=v(nrn_monit);
	if savevideo & mod(t,spiking_bin) == 0
		myMovie = heatmap(ncells, gc_firings, t, skip_t, h, myMovie, ccol, spiking_bin);
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

disp(size(find(gc_firings(:,2)==nrn_monit))); % number of spikes from selected neuron
%disp("voltage");
%disp(voltage);
if false
	figure(1);
	plot(gc_voltage);
	ylim([-80 30])
	title("gc voltage")
	figure(2);
	plot(in_voltage);
	ylim([-80 30])
	title("in voltage")
end

function [v, u] = iznrn(v, u, p, fired, Ie, Ii)
	a=p(:,1);b=p(:,2);c=p(:,3);d=p(:,4);
	v(fired)=c(fired);
	u(fired)=u(fired)+d(fired);
	v=v+(0.04*v.^2+5*v+140-u+Ie-Ii); % step 1.0 ms
	u=u+a.*(b.*v-u);

	% avoid NaN issue
	nans = find(isnan(v));
	if size(nans)>0
		v(nans) = c(1);
	end
	% hard limit on v to avoid inf value issue
	highval = find(v>1000000);
	if size(highval)>0
		v(highval) = c(1);
	end
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
	% report spike times in a bin of time. detect current spikes in current ms.
	spike_found = false;
	all_spike_times = (find(firings(:,2)==ni));
	for si=1:size(all_spike_times)
		spike_time = firings(all_spike_times(si),1);
        if spike_time == t
            spike_found = true;
        end
    end
end

function [gc_ie, vi, ui] = gc_in_signal(gc_ie, t, gc_firings, in_fired, vi, ui, p2, gcintau, ncells, nrn_monit, gc_to_in_wt)
    % generate gc to in signaling
    global cap_ue tau_ue tau_xe tau_de gei u_ei x_ei;
	gc_firing = zeros(ncells,1);
	for i=1:ncells
        spike_found = find_spike(i,t,gc_firings);
		if spike_found == true
	        gc_firing(i) = gc_firing(i)+1;
	    end    
    end
	% simple synapse for one-to-one connections
	if 0 
		weight = gc_to_in_wt*gc_firing;
	    gc_ie = gc_ie + -gc_ie/gcintau + weight/gcintau; %gc_ie = gc_ie + (-gc_ie + weight)/gcintau;
	end
    % tm model synapse
    if 1
    	%gc_to_in_wt = 100;
    	weights = gc_to_in_wt*gc_firing;
    	%[u_ei x_ei gc_ie] = tm_synapse(u_ei,x_ei,gc_ie,cap_ue,tau_ue,tau_xe, ...
    	%					tau_de,gei,gc_to_in_wt,gc_firing);
    	[u_ei x_ei gc_ie] = tm_synapse(u_ei,x_ei,gc_ie,cap_ue,tau_ue,tau_xe, ...
    						tau_de,gei,weights,gc_firing);
    	%ts=0.5; % time step of 0.5 ms for numerical stability
		%u_ei=u_ei+ts*((-u_ei/tau_u)+(cap_u*(1-u_ei))*spk);
		%disp("u_ei");
		%disp(size(u_ei));
		%disp(size(spk));
		%u_ei=u_ei+ts*((-u_ei/tau_u)+(cap_u*(1-u_ei)).*spk');
		%disp(size(test));
    end

    %disp(gc_firing(772));
   	%disp(gc_ie(772));
   	%fprintf("t:%d i:%f\n",t,gc_ie(772));
	gc_ii = zeros(1,900)';
	%n= ncells;
	%fprintf("t:%d %f=%f+(0.04*%f.^2+5*%f+140-%f+%f-%f)\n",t,(vi(n)+(0.04*vi(n).^2+5*vi(n)+140-ui(n)+gc_ie(n)-gc_ii(n))),vi(n),vi(n),vi(n),ui(n),gc_ie(n),gc_ii(n));	
	[vi, ui] = iznrn(vi, ui, p2, in_fired, gc_ie, gc_ii);
	%fprintf("t:%d vi:%f\n",t,vi(nrn_monit));
end

function [in_ii, in_firings] = in_gc_signal(t, mex_hat, in_firings, ncells, in_ii, gcintau, nrn_monit, in_to_gc_wt);
	% generate in to gc signaling
	global cap_ui tau_ui tau_xi tau_di gie u_ie x_ie;
	in_firing = zeros(ncells,1);
	for i=1:ncells
        spike_found = find_spike(i,t,in_firings);
		if spike_found == true
	        in_firing(i) = in_firing(i)+1;
	    end    
    end
    if 0 % simple synapse for one-to-many connections
    	% 900x900 matrix converted to 900x1 matrix because it appears calcs are 
    	% equivalent to larger matrix and saves comp time. E.g., all individual 
    	% weights are added up for 30x30 gc layer eventually anyway.
		ingc_current = ((mex_hat*in_to_gc_wt).*in_firing');
		ingc_summed = ingc_current*ones(ncells,1);  
	    in_ii = in_ii - in_ii/gcintau + ingc_summed/gcintau;		
	end
	if 1 % tm synapse model
		%weights = mex_hat*in_to_gc_wt;
		weights = ((mex_hat*in_to_gc_wt).*in_firing');
		weights = weights*ones(ncells,1);
		%in_ii = in_ii - in_ii/30 + weights.*0.033;%.*in_firing;
		[u_ie x_ie in_ii] = tm_synapse(u_ie,x_ie,in_ii,cap_ui,tau_ui,tau_xi, ...
    						tau_di,gie,weights,in_firing);
	end
	%fprintf("t:%d i:%d\n",t,sum(find(in_firing(1,:)>=1),1));
	%fprintf("t:%d s:%d\n",t,sum(in_firing(1,:)));
	%disp(in_firing(1,:));
    %fprintf("t:%d i:%f\n",t,in_ii(nrn_monit));
end

function Ii = inhib_curr(Ii, t, mex_hat, firings, tau, nrn_monit)
	% generate inhibitory currents
	gc_firing = zeros(size(mex_hat,1)); 
	for i=1:size(Ii)
        spike_found = find_spike(i,t,firings);
		if spike_found == true
	        gc_firing(:,i) = gc_firing(:,i)+1;
	    end    
    end
    in_current = (mex_hat*gc_firing')';
    fprintf("t:%d s:%d\n",t,sum(gc_firing(1,:)));

	% calculate tau factor
	o = ones(size(mex_hat(:,1)));
	in_summed = in_current'*o; in_summed2 = in_summed;
	in_summed2 = in_summed2.*(in_summed2>0); % no negative values
	Ii = Ii + (in_summed2 - Ii)/tau;
	%fprintf("t:%d i:%f\n",t,Ii(nrn_monit));
end

function myMovie = heatmap(ncells, firings, t, skip_t, h, myMovie, ccol, spiking_bin)
	binned_firing = [];
	for i=1:sqrt(ncells)
		temp = [];
		for j=1:sqrt(ncells)
			nrn_i = ((i-1)*sqrt(ncells))+j;
			spk_t = fbin(nrn_i, t-spiking_bin, t, firings);
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
    %caxis([0 4.0])
    caxis([0 8.0])
	cb = colorbar;
	%set(cb, 'ylim', [0 5.5]); % set colorbar range
	set(cb, 'ylim', [0 8.5]); % set colorbar range
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