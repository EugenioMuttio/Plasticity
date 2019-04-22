function [sigma_n1,int_vars_n1,dev_sigma_n1] = maps_plasJ2(matprop,Ce,sigma_n,eps_rate,int_vars_n,delta_t)

%Material Properties
E=matprop(1);
nu = matprop(9);
mu = E/(2*(1+nu));
sigma_y=matprop(2);
hard_type=matprop(3);
K=matprop(4);
H=matprop(5);
delta=matprop(6);

%Internal Variables from last step
int_vars_n1=int_vars_n;
eps_n=int_vars_n(:,1);
eps_n1=int_vars_n(:,2);
eps_p_n=int_vars_n(:,3);
eps_p_n1=int_vars_n(:,4);
xi_n=int_vars_n(:,5);
xi_n1=int_vars_n(:,6);
xibar_n=int_vars_n(:,7);
xibar_n1=int_vars_n(:,8);
q_n=int_vars_n(:,11);
qbar_n=int_vars_n(:,13);
q_n1=q_n;
qbar_n1=qbar_n;

%Elastic Sigma using the GIVEN strains
sigma_nt=Ce*eps_n;
sigma_n1t=Ce*eps_n1;

%Sigma trial is considered as the elastic stress by using the GIVEN STRAINS
% at time (n+1) minus the plastic strain COMPUTED at the step before
% which initially is zero
sigma_trial_n=Ce*(eps_n-eps_p_n);
sigma_trial=Ce*(eps_n1-eps_p_n);
q_trial=q_n(1);%-K*xi_n(1);%scalar
qbar_trial = -2/3*H*xibar_n;%vector

shp = (1/3)*(sigma_trial(1)+sigma_trial(2)+sigma_trial(3));
dev_sigma_trial_n = sigma_trial_n-[shp shp shp 0 0 0]';
dev_sigma_trial = sigma_trial-[shp shp shp 0 0 0]';

n_trial = ((dev_sigma_trial-qbar_trial)/norm(dev_sigma_trial-qbar_trial));

%ELASTIC Sigma rate used only for the sign, to know when the rate 
%changes from load -> unload
sigma_rate_trial=(sigma_n1t-sigma_nt)/delta_t;
dev_sigma_rate_trial=(dev_sigma_trial-dev_sigma_trial_n)/delta_t;

%Parameter to include in the plasticity models
%gamma_n=int_vars_n(9); 
gamma_n1=int_vars_n(10);

if hard_type==0
    H=0;
    K=0;
    hard_type=1;
end

if hard_type==1
    
%-- LINEAR ISOTROPIC / KINEMATIC HARDENING -------------------------------
%-------------------------------------------------------------------------
    sigma_lim=sigma_y-q_trial;
    norm_sigm = norm(dev_sigma_trial-qbar_trial);
    ftrial = norm_sigm-sqrt(2/3)*sigma_lim;
    % Sigma condition f(sigma)=0
    if abs(round(norm_sigm,8))<=(round(sqrt(2/3)*sigma_lim,8))  
        sigma_n1=sigma_trial;
        eps_p_n1=eps_p_n;
        
        %Internal variables
        xi_n1=xi_n;
        xibar_n1=xibar_n;
        
        q_n1=q_n;
        qbar_n1=qbar_n;
    else
        % Sigma rate condition df/d(sigma)*sigma_trial>0
        if (n_trial'*dev_sigma_rate_trial)>0
            %Plastic Multiplier
            gamma_n1=((2*mu+2/3*K+2/3*H)^(-1))*ftrial;
        else
            gamma_n1=0;
        end
        
        sigma_n1 = sigma_trial-gamma_n1*2*mu*n_trial;
        %sigma_n1 = sigma_trial-gamma_n1*Ce*n_trial;
        q_n1(1) = q_trial-gamma_n1*K*sqrt(2/3);
        qbar_n1 = qbar_trial+gamma_n1*2/3*H*n_trial;
        
        %Internal variables (n+1) computation
        eps_p_n1 = eps_p_n+gamma_n1*n_trial;
        xi_n1 = [xi_n(1)+gamma_n1*sqrt(2/3) zeros(1,5)]';
        xibar_n1 = xibar_n-gamma_n1*n_trial;
    end
    
    
elseif hard_type==2
%-- EXPONENTIAL SATURATION LAW + LINEAR HARDENING ------------------------
%-------------------------------------------------------------------------
    sigma_inf=1.5*sigma_y;
    %Exponential Saturation Law Function for NR
    func=@(xi,sigma_inf,sigma_y,delta,K)(sigma_inf-sigma_y)*(1-exp(-delta*xi))+K*xi;
    %Derivative Exponential Saturation Law Function for NR
    dfdxi=@(xi,delta,K,sigma_inf,sigma_y)(sigma_inf-sigma_y)*(delta*exp(-delta*xi))+K;
    
    %Exponential saturation law
    sigma_lim=sigma_y-q_trial;
    norm_sigm = norm(dev_sigma_trial-qbar_trial);
    ftrial = norm_sigm-sqrt(2/3)*sigma_lim;
    
    % Sigma condition f(sigma)=0
    if abs(round(norm_sigm,8))<=(round(sqrt(2/3)*sigma_lim,8))  
        sigma_n1=sigma_trial;
        eps_p_n1=eps_p_n;
        
        %Internal variables
        xi_n1=xi_n;
        xibar_n1=xibar_n;
        q_n1=q_n;
        qbar_n1=qbar_n;
    else
         % Sigma rate condition df/d(sigma)*sigma_trial>0
        if (n_trial'*dev_sigma_rate_trial)>0
            %Material parameters for evaluation inside NR
            mat_param=[xi_n(1) mu H K delta sigma_y sigma_inf];
            %Newton-Rapshon
            gamma_n1=NR_J2(func,dfdxi,ftrial,mat_param);
            
        else
            gamma_n1=0;            
        end
        
        
        %Internal variables (n+1) computation
        xi_n1(1)=xi_n(1)+sqrt(2/3)*gamma_n1;
        xibar_n1=xibar_n+gamma_n1*sign(sigma_trial-qbar_n);
        q_n1(1)=q_trial-func(xi_n1(1),sigma_inf,sigma_y,delta,K)+func(xi_n(1),sigma_inf,sigma_y,delta,K);
        qbar_n1=qbar_trial+gamma_n1*2/3*H*n_trial;
        
        sigma_n1=sigma_trial-gamma_n1*2*mu*n_trial;
        eps_p_n1=eps_p_n+gamma_n1*n_trial;
        
    end
    
    
end
shp_n1 = (1/3)*(sigma_n1(1)+sigma_n1(2)+sigma_n1(3));
dev_sigma_n1 = sigma_n1-[shp_n1 shp_n1 shp_n1 0 0 0]';


%%UPDATE VALUES%%
int_vars_n1(:,1)=eps_n;
int_vars_n1(:,2)=eps_n1;
int_vars_n1(:,3)=eps_p_n;
int_vars_n1(:,4)=eps_p_n1;
int_vars_n1(:,5)=xi_n;
int_vars_n1(:,6)=xi_n1;
int_vars_n1(:,7)=xibar_n;
int_vars_n1(:,8)=xibar_n1;
int_vars_n1(1,10)=gamma_n1;
int_vars_n1(:,12)=q_n1;
int_vars_n1(:,14)=qbar_n1;











