cap_u = 0.2; % U, utilization
tau_u = 20.0; % U signal decay time constant
tau_x = 100.0; % x signal decay time constant
tau_d = 5.0; % x signal decay time constant
g = 1.0;
u = 1; % u before spike update
x = 1; % x before spike update
i = 1;
A = 2.0;
spk = 1; % presence of a spike triggering the synapse. 0=no spk, 1=spk present
%tspan = [0 1]; % time span to compute (1 ms)

%[t_u,du_dt] = ode23(@(t,u) u_deriv(t,u,tau_u,cap_u,spk), tspan, u_pre);
%u_post = du_dt(end);
%[t_x,dx_dt] = ode23(@(t,x) x_deriv(t,x,tau_x,u_post,spk), tspan, x_pre);
%x_post = dx_dt(end);
%[t_i,di_dt] = ode23(@(t,i) i_deriv(t,i,tau_d,A,u_post,x_post,spk), tspan, i_pre);
%i_post = di_dt(end)

tm_synapse(u,x,i,cap_u,tau_u,tau_x,tau_d,g,A,spk)