function [sigma_n1,int_vars_n1] = maps_plas(matprop,sigma_n,eps_rate,int_vars_n,delta_t)

%Material Properties
E=matprop(1);
sigma_y=matprop(2);
hard_type=matprop(3);
K=matprop(4);
H=matprop(5);
delta=matprop(6);
eta=matprop(8);

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

%Sigma trial is considered as the elastic stress computed by using
%the GIVEN strains
sigma_trial=sigma_n1t;

%ELASTIC Sigma rate used only for the sign, to know when the rate 
%changes from load -> unload
sigma_rate_trial=(sigma_n1t-sigma_nt)/delta_t;

%Parameter to include in the plasticity models
gamma_n=int_vars_n(9); 
gamma_n1=0;

if hard_type==0
%-- Perfect Plasticity ---------------------------------------------------
%-------------------------------------------------------------------------
    %Plastic strain computation - Backward Euler (BE) time integration
    eps_p_rate=gamma_n1*sign(sigma_trial);
    eps_p_n1=eps_p_n+eps_p_rate*delta_t;

    %Return Mapping Algorithm
    %Sigma rate computed as Saracibar slide 1 page 30 
    sigma_rate=sigma_rate_trial-E*eps_p_rate;

    %New Sigma computed by considering the plastic strains
    sigma_n1=(sigma_rate)*delta_t+sigma_n;

    % Sigma condition f(sigma)=0
    % sigma trial is the previous one (**maybe this can change to optimize
    % the code implemented by using the sigma of this step)
    if abs(round(sigma_n1,8))>=sigma_y
        % Sigma rate condition df/d(sigma)*sigma_trial>0
        if (sigma_rate_trial)*sign(sigma_trial)>0
            gamma_n1=eps_rate*sign(sigma_trial);
        else
            gamma_n1=0;
        end
        
        sigma_n1=sigma_y*sign(sigma_trial);
    end

elseif hard_type==1
%-- LINEAR ISOTROPIC / KINEMATIC HARDENING -------------------------------
%-------------------------------------------------------------------------
    sigma_trial=E*(eps_n1-eps_p_n);
    q=-K*xi_n;
    qbar=H*xibar_n;
    sigma_lim=sigma_y-q;
    
    % Sigma condition f(sigma)=0
    if abs(round(sigma_trial-qbar,8))<=(round(sigma_lim,8))  
        sigma_n1=sigma_trial;
        eps_p_n1=eps_p_n;
    else
        ftrial=abs(sigma_trial-qbar)-sigma_y+q;
        if (sigma_rate_trial)*sign(sigma_trial)>0
            %Plastic Multiplier
            gamma_n1=((1/delta_t)*(E+K+H+eta/delta_t)^(-1))*ftrial;
        else
            gamma_n1=0;
        end
        
        sigma_n1=sigma_trial-gamma_n1*delta_t*E*sign(sigma_trial-qbar);
        eps_p_n1=eps_p_n+gamma_n1*delta_t*sign(sigma_trial-qbar);
         
    end

    %Internal variables (n+1) computation
    xi_rate=abs(gamma_n1);
    xi_n1=xi_n+xi_rate*delta_t;
    xibar_rate=gamma_n1*sign(sigma_n-qbar);
    xibar_n1=xibar_n+xibar_rate*delta_t;
    
elseif hard_type==2
%-- EXPONENTIAL SATURATION LAW + LINEAR HARDENING ------------------------
%-------------------------------------------------------------------------
    sigma_inf=200*sigma_y;
    %Exponential Saturation Law Function for NR
    func=@(xi,sigma_inf,sigma_y,delta,K)(sigma_inf-sigma_y)*(1-exp(-delta*xi))+K*xi;
    %Derivative Exponential Saturation Law Function for NR
    dfdxi=@(xi,delta,K)(delta*exp(-delta*xi))+K;
    
    %Exponential saturation law
    q=-func(xi_n,sigma_inf,sigma_y,delta,K);
    qbar=H*xibar_n;
    sigma_lim=sigma_y-q;
    
    % Sigma condition f(sigma)=0
    if abs(round(sigma_n-qbar,8))>=(round(sigma_lim,8))
         
        sigma_n=sigma_y*sign(sigma_trial-qbar)-q*sign(sigma_trial-qbar)+qbar;
        % Sigma rate condition df/d(sigma)*sigma_trial>0
        if (sigma_rate_trial)*sign(sigma_trial-qbar)>0
            %Nonlinear Problem, then Newton-Raphson method is employed to
            %compute the plastic multiplier gamma 
            
            ftrial=(abs(sigma_trial-qbar)-sigma_y+q)*sign(sigma_n-qbar);
            %Material parameters for evaluation inside NR
            mat_param=[xi_n E H K delta sigma_y sigma_inf];
            %Newton-Rapshon
            gamma_n1=NR_1D(func,dfdxi,ftrial,mat_param);
            
        else
            gamma_n1=0;
        end
        
    end
    
    
    %Plastic strain computation - Backward Euler (BE) time integration
    eps_p_rate=gamma_n1*sign(sigma_n);
    eps_p_n1=eps_p_n+eps_p_rate*delta_t;

    %Return Mapping Algorithm
    %Sigma rate computed as Saracibar slide 1 page 30 
    sigma_rate=sigma_rate_trial-E*eps_p_rate;

    %New Sigma computed by considering the plastic strains
    sigma_n1=(sigma_rate)*delta_t+sigma_n;
    
    %Internal variables (n+1) computation
    xi_rate=abs(gamma_n1);
    xi_n1=xi_n+xi_rate*delta_t;
    xibar_rate=gamma_n1*sign(sigma_n1-qbar);
    xibar_n1=xibar_n+xibar_rate*delta_t;
        
    %Exponential saturation law
    q=-func(xi_n1,sigma_inf,sigma_y,delta,K);
    qbar=H*xibar_n1;
    sigma_lim=sigma_y-q;
    
    % Sigma condition f(sigma)=0
    if abs(round(sigma_n1-qbar,8))>=(round(sigma_lim,8))
        sigma_n1=sigma_y*sign(sigma_trial-qbar)-q*sign(sigma_trial-qbar)+qbar;
    end
    
end



%%UPDATE VALUES%%
int_vars_n1(1)=eps_n;
int_vars_n1(2)=eps_n1;
int_vars_n1(3)=eps_p_n;
int_vars_n1(4)=eps_p_n1;
int_vars_n1(5)=xi_n;
int_vars_n1(6)=xi_n1;
int_vars_n1(7)=xibar_n;
int_vars_n1(8)=xibar_n1;
int_vars_n1(10)=gamma_n1;










