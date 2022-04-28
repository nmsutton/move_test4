function [u x i]=tm_synapse(u,x,i,cap_u,tau_u,tau_x,tau_d,g,A,spk)
	% perform numerical integration of the tm synapse model for 1 ms
	% this is the CARLsim equations version
	ts=1;%0.5; % time step of 0.5 ms for numerical stability
	u=u+ts*((-u/tau_u)+(cap_u*(1-u)).*spk);
	%u=u+ts*((-u/tau_u)+(cap_u*(1-u)).*spk);
	x=x+ts*(((1-x)/tau_x)-u.*x.*spk);
	%x=x+ts*(((1-x)/tau_x)-u.*x.*spk);
	%disp("ux");
	%disp(size(u));
	%disp(size(x));
	%u=ones(900,1);%*0.00667;
	%x=ones(900,1)*.033;%.021;
	%i=i+ts*(-i/(tau_d*g)+A*u.*x.*spk);
	i=i+ts*(-i/(tau_d*g)+A*u.*x.*spk);