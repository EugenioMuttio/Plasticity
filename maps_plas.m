function [sigma_n1,int_vars_n1] = maps_plas(matprop,eps_rate,int_vars_n,delta_t)

E=matprop(1);
sigma_y=matprop(2);
hard_type=matprop(3);
K=matprop(4);

int_vars_n1=int_vars_n;
eps_n=int_vars_n(1);
eps_n1=int_vars_n(2);
eps_p_n=int_vars_n(3);
eps_p_n1=int_vars_n(4);
xi_n=int_vars_n(5);
xi_n1=int_vars_n(6);
xibar_n=int_vars_n(7);
xibar_n1=int_vars_n(8);

sigma_n=E*eps_n;
sigma_n1=E*eps_n1;
sigma_rate=(sigma_n1-sigma_n)/delta_t;

%Sigma GIVEN from the user to compare and test
sigma_trial=sigma_n1;
%Parameter to include in the plasticity models
gamma=0; 
%Sigma COMPUTED from the plasticity model
sigma=sigma_n1;



if hard_type==0
%-- Perfect Plasticity ---------------------------------------------------
%-------------------------------------------------------------------------
    if abs(sigma_trial)>sigma_y
        gamma=eps_rate*sign(sigma_rate);
        sigma=sigma_y*sign(sigma_trial);
    else
        gamma=0;

    end

    eps_p_rate=gamma*sign(sigma_trial);
    eps_p_n1=eps_p_rate*delta_t+eps_p_n;

elseif hard_type==1
%-- LINEAR HARDENING -----------------------------------------------------
%-------------------------------------------------------------------------
    q=-K*xi_n;
    sigma_lim=sigma_y-q;
    if abs(sigma_trial)>(sigma_lim)
        sigma_trial=sigma_y*sign(sigma_trial)-q*sign(sigma_trial);
        gamma=E*eps_rate*sign(sigma_rate)/(E+K);
        sigma=(sigma_rate-gamma*sign(sigma_trial))*delta_t+sigma_trial;
        
    else
        gamma=0;

    end
    
    eps_p_rate=gamma*sign(sigma_trial);
    eps_p_n1=eps_p_rate*delta_t+eps_p_n;
    
    xi_rate=abs(gamma);
    xi_n1=xi_rate*delta_t+xi_n;
    
    
end

%%UPDATE VALUES%%
sigma_n1=sigma;
int_vars_n1(1)=eps_n;
int_vars_n1(2)=eps_n1;
int_vars_n1(3)=eps_p_n;
int_vars_n1(4)=eps_p_n1;
int_vars_n1(5)=xi_n;
int_vars_n1(6)=xi_n1;
int_vars_n1(7)=xibar_n;
int_vars_n1(8)=xibar_n1;










