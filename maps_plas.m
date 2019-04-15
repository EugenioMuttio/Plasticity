function [eps_p_n1,sigma_n1] = maps_plas(matprop,eps_rate,eps_n,eps_n1,eps_p_n,eps_p_n1,delta_t)

E=matprop(1);
sigma_y=matprop(2);

sigma_n=E*eps_n;
sigma_n1=E*eps_n1;
sigma_rate=(sigma_n1-sigma_n)/delta_t;

sigma_rate_trial=E*eps_rate;
gamma=0; 
sigma=sigma_n1;
f=abs(sigma)-sigma_y;
eps_p_rate=0;
if abs(sigma_n1)>sigma_y
    sigma=sigma_y*sign(sigma_n1);
    gamma=eps_rate*sign(sigma_rate);
else
    gamma=0;

end

eps_p_rate=gamma*sign(sigma_n1);
eps_p_n1=eps_p_rate*delta_t+eps_p_n;

sigma_n1=sigma;










