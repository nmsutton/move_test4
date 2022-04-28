function i_post=tm_synapse()
	[t_u,du_dt] = ode23(@(t,u) u_deriv(t,u,tau_u,cap_u,spk), tspan, u_pre);
	u_post = du_dt(end);
	[t_x,dx_dt] = ode23(@(t,x) x_deriv(t,x,tau_x,u_post,spk), tspan, x_pre);
	x_post = dx_dt(end);
	[t_i,di_dt] = ode23(@(t,i) i_deriv(t,i,tau_d,A,u_post,x_post,spk), tspan, i_pre);
	i_post = di_dt(end);