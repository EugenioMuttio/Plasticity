function [sigma_n1,int_vars_n1] = maps_plas(matprop,sigma_n,eps_rate,int_vars_n,delta_t)

%Material Properties
E=matprop(1);
sigma_y=matprop(2);
hard_type=matprop(3);
K=matprop(4);
H=matprop(5);

%Internal Variables from last step
int_vars_n1=int_vars_n;
eps_n=int_vars_n(1);
eps_n1=int_vars_n(2);
eps_p_n=int_vars_n(3);
eps_p_n1=int_vars_n(4);
xi_n=int_vars_n(5);
xi_n1=int_vars_n(6);
xibar_n=int_vars_n(7);
xibar_n1=int_vars_n(8);

%Elastic Sigma using the GIVEN strains
sigma_nt=E*eps_n;
sigma_n1t=E*eps_n1;

%Elastic Sigma rate used only for the sign, to know when the rate 
%changes from load -> unload
sigma_rate_trial=(sigma_n1t-sigma_nt)/delta_t;


%Sigma trial is the corresponding to the PREVIOUS step
%Starts with zero, then changes to elastic and then depends on the model.
sigma_trial=sigma_n;
%Parameter to include in the plasticity models
gamma=0; 


if hard_type==0
%-- Perfect Plasticity ---------------------------------------------------
%-------------------------------------------------------------------------
    % Sigma condition f(sigma)=0
    % sigma trial is the previous one (**maybe this can change to optimize
    % the code implemented by using the sigma of this step)
    if abs(sigma_trial)>=sigma_y
        sigma_n=sigma_y*sign(sigma_trial);
    
        % Sigma rate condition df/d(sigma)*sigma_trial>0
        if (sigma_rate_trial)*sign(sigma_n)>0
            gamma=eps_rate*sign(sigma_trial);
        else
            gamma=0;
        end
    end
    
    %Plastic sttrain computation
    eps_p_rate=gamma*sign(sigma_trial);
    eps_p_n1=eps_p_rate*delta_t+eps_p_n;
    
    %Sigma rate computed as Saracibar slide 1 page 30 
    sigma_rate=sigma_rate_trial-E*eps_p_rate;
    
    %New Sigma computed by considering the plastic strains
    sigma_n1=(sigma_rate)*delta_t+sigma_n;
    
    %Sigma model verification (**this part could be modified to optimize)
    if abs(sigma_n1)>sigma_y
        sigma_n1=sigma_y*sign(sigma_trial);
    end

elseif hard_type==1
%-- LINEAR HARDENING -----------------------------------------------------
%-------------------------------------------------------------------------
    q=-K*xi_n;
    qbar=H*xibar_n;
    sigma_lim=sigma_y-q;
    % Sigma condition f(sigma)=0
    if abs(round(sigma_trial-qbar,8))>=(round(sigma_lim,8))
        sigma_n=sigma_y*sign(sigma_trial-qbar)-q*sign(sigma_trial-qbar)+qbar;
        
        % Sigma condition f(sigma)=0
        if (sigma_rate_trial)*sign(sigma_n-qbar)>0
            gamma=E*eps_rate*sign(sigma_trial-qbar)/(E+K+H);
        else
            gamma=0;
        end
    end
    
    %Plastic sttrain computation
    eps_p_rate=gamma*sign(sigma_trial-qbar);
    eps_p_n1=eps_p_rate*delta_t+eps_p_n;

    %Sigma rate computed as Saracibar slide 1 page 30 
    sigma_rate=sigma_rate_trial-E*eps_p_rate;

    %New Sigma computed by considering the plastic strains
    sigma_n1=(sigma_rate)*delta_t+sigma_n;

    %Internal variables computation
    xi_rate=abs(gamma);
    xi_n1=xi_rate*delta_t+xi_n;
    xibar_rate=gamma*sign(sigma_n-qbar);
    xibar_n1=xibar_rate*delta_t+xibar_n;
    
    %Sigma model verification (**this part could be modified to optimize)
    q=-K*xi_n1;
    qbar=H*xibar_n1;
    sigma_lim=sigma_y-q;
    
    if abs(sigma_n1-qbar)>sigma_lim
        sigma_n1=sigma_y*sign(sigma_trial-qbar)-q*sign(sigma_trial-qbar)+qbar;
    end
    
    
end



%%UPDATE VALUES%%
%sigma_n1=sigma;
int_vars_n1(1)=eps_n;
int_vars_n1(2)=eps_n1;
int_vars_n1(3)=eps_p_n;
int_vars_n1(4)=eps_p_n1;
int_vars_n1(5)=xi_n;
int_vars_n1(6)=xi_n1;
int_vars_n1(7)=xibar_n;
int_vars_n1(8)=xibar_n1;










