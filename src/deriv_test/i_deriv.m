function di_dt=i_deriv(t,i,tau_d,A,u_post,x_post,spk)
	di_dt=(-i/tau_d)+A*u_post*x_post*spk;