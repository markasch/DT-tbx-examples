function dydt = ode_unstable(t,y)
dydt = -(1/t.^2) + 10.0*(y-1./t);