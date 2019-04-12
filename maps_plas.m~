function sigma_n1 = maps_plas(matprop,eps_rate,eps_n,eps_n1)

E=matprop(1);
sigma_y=matprop(2);

sigma_n=E*eps_n;
sigma_trial=E*eps_n1;



%Perfect Plasticity
if abs(sigma_trial)>sigma_y
    eps_p=eps_n1-sigma_y/E;
    sigma_trial=sigma_y*sign(sigma_trial);
    gamma=abs(eps_rate);
else
    gamma=0;
end

sigma_n1=sigma_trial;





