function dx_dt=x_deriv(t,x,tau_x,u_post,spk)
	dx_dt=((1-x)/tau_x)-u_post*x*spk;