function du_dt=u_deriv(t,u,tau_u,cap_u,spk)
	du_dt=(-u/tau_u)+(cap_u*(1-u))*spk;