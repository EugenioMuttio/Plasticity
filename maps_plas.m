function [sigma_n1,int_vars_n1] = maps_plas(matprop,sigma_n,eps_rate,int_vars_n,delta_t)

%Material Properties
E=matprop(1);
sigma_y=matprop(2);
hard_type=matprop(3);
K=matprop(4);
H=matprop(5);
delta=matprop(6);

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
q_n=int_vars_n(11);
qbar_n=int_vars_n(13);
q_n1=q_n;
qbar_n1=qbar_n;

%Elastic Sigma using the GIVEN strains
sigma_nt=E*eps_n;
sigma_n1t=E*eps_n1;

%Sigma trial is considered as the elastic stress by using the GIVEN STRAINS
% at time (n+1) minus the plastic strain COMPUTED at the step before
% which initially is zero
sigma_trial=E*(eps_n1-eps_p_n);

%ELASTIC Sigma rate used only for the sign, to know when the rate 
%changes from load -> unload
sigma_rate_trial=(sigma_n1t-sigma_nt)/delta_t;

%Parameter to include in the plasticity models
gamma_n=int_vars_n(9); 
gamma_n1=0;

if hard_type==0
%-- Perfect Plasticity ---------------------------------------------------
%-------------------------------------------------------------------------

    % Sigma condition f(sigma)=0
    if abs(round(sigma_trial,8))<=sigma_y
        sigma_n1=sigma_trial;
        eps_p_n1=eps_p_n;
    else        
        ftrial=abs(sigma_trial)-sigma_y;
        % Sigma rate condition df/d(sigma)*sigma_trial>0
        if (sigma_rate_trial)*sign(sigma_trial)>0
            %Plastic Multiplier
            gamma_n1=((E)^(-1))*ftrial;
        else
            gamma_n1=0;
        end
        
        sigma_n1=sigma_trial-gamma_n1*E*sign(sigma_trial);
        eps_p_n1=eps_p_n+gamma_n1*sign(sigma_trial);
         
    end

elseif hard_type==1
%-- LINEAR ISOTROPIC / KINEMATIC HARDENING -------------------------------
%-------------------------------------------------------------------------
    q=-K*xi_n;
    qbar=H*xibar_n;
    sigma_lim=sigma_y-q;
    % Sigma condition f(sigma)=0
    if abs(round(sigma_trial-qbar,8))<=(round(sigma_lim,8))  
        sigma_n1=sigma_trial;
        eps_p_n1=eps_p_n;
        
        %Internal variables
        xi_n1=xi_n;
        xibar_n1=xibar_n;
        
        q_n1=q_n;
        qbar_n1=qbar_n;
    else
        ftrial=abs(sigma_trial-qbar)-sigma_y+q;
        % Sigma rate condition df/d(sigma)*sigma_trial>0
        if (sigma_rate_trial)*sign(sigma_trial-qbar)>0
            %Plastic Multiplier
            gamma_n1=((E+K+H)^(-1))*ftrial;
        else
            gamma_n1=0;
        end
        
        sigma_n1=sigma_trial-gamma_n1*E*sign(sigma_trial-qbar);
        eps_p_n1=eps_p_n+gamma_n1*sign(sigma_trial-qbar);
        
        %Internal variables (n+1) computation
        xi_n1=xi_n+gamma_n1;
        xibar_n1=xibar_n+gamma_n1*sign(sigma_trial-qbar);
        
        qbar_n1=qbar_n+gamma_n1*H*sign(sigma_trial-qbar_n);
        q_n1=q_n-gamma_n1*K;
    end
    


elseif hard_type==2
%-- EXPONENTIAL SATURATION LAW + LINEAR HARDENING ------------------------
%-------------------------------------------------------------------------
    %if eps_rate*sign(sigma_trial)>0
        sigma_inf=1.2*sigma_y;
    %else
     %   sigma_inf=1.5*sigma_y;
    %end
    %Exponential Saturation Law Function for NR
    func=@(xi,sigma_inf,sigma_y,delta,K)(sigma_inf-sigma_y)*(1-exp(-delta*xi))+K*xi;
    %Derivative Exponential Saturation Law Function for NR
    dfdxi=@(xi,delta,K,sigma_inf,sigma_y)(sigma_inf-sigma_y)*(delta*exp(-delta*xi))+K;
    
    %Exponential saturation law
    q=-func(xi_n,sigma_inf,sigma_y,delta,K);
    qbar=H*xibar_n;
    sigma_lim=sigma_y-q_n;
    
    % Sigma condition f(sigma)=0
    if abs(round(sigma_trial-qbar_n,8))<=(round(sigma_lim,8))  
        sigma_n1=sigma_trial;
        eps_p_n1=eps_p_n;
        
        %Internal variables
        xi_n1=xi_n;
        xibar_n1=xibar_n;
        q_n1=q_n;
        qbar_n1=qbar_n;
    else
        ftrial=abs(sigma_trial-qbar_n)-sigma_y+q_n;
        % Sigma rate condition df/d(sigma)*sigma_trial>0
        if (sigma_rate_trial)*sign(sigma_trial-qbar)>0
            %Material parameters for evaluation inside NR
            mat_param=[xi_n E H K delta sigma_y sigma_inf];
            %Newton-Rapshon
            gamma_n1=NR_1D(func,dfdxi,ftrial,mat_param);
            
        else
            gamma_n1=0;            
        end
        
        
        %Internal variables (n+1) computation
        xi_n1=xi_n+gamma_n1;
        xibar_n1=xibar_n+gamma_n1*sign(sigma_trial-qbar_n);
        q_n1=q_n-func(xi_n1,sigma_inf,sigma_y,delta,K)+func(xi_n,sigma_inf,sigma_y,delta,K);
        qbar_n1=qbar_n+gamma_n1*H*sign(sigma_trial-qbar_n);
        
        sigma_n1=sigma_trial-gamma_n1*E*sign(sigma_trial-qbar_n);
        eps_p_n1=eps_p_n+gamma_n1*sign(sigma_trial-qbar_n);
        
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
int_vars_n1(12)=q_n1;
int_vars_n1(14)=qbar_n1;











