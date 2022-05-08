function [u x A i]=tm_synapse(u,x,A,i,cap_u,tau_u,tau_x,tau_d,g,w,spk)
	% perform numerical integration of the tm synapse model for 1 ms
	% this is the simplified (Tsodyks, Pawelzik, and Markram, 1998) 
	% equations version from (Moradi et al., 2022).

	
    u=u+(-u/tau_u)+(cap_u*(1-u)).*spk;
    x=x+((1-x-A)/tau_x)-u.*x.*spk;
	%u=ones(900,1);%1;%ones(900,1);%*0.00667;
	%x=ones(900,1)*.0333;%.0333;%ones(900,1)*.0333;%.021;    
    A=A+(-A/tau_d)+u.*x.*spk;
    %A=ones(900,1)*.0333;%.0333;%ones(900,1)*.0333;%.021;
    i=g.*w.*A;
    
    %{
    ts=0.025;
    i = 0;
    for i=1:(1/ts)       
        u=u+ts*((-u/tau_u)+(cap_u*(1-u)).*spk);
        x=x+ts*(((1-x-A)/tau_x)-u.*x.*spk);
        A=A+ts*((-A/tau_d)+u.*x.*spk);
        i=i+ts*(g.*w.*A);
    end
    %}