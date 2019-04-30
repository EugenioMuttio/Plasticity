clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program for modelling plasticity damage model
% (Elemental gauss point level)
% -------------------------------
% Developed by:
% Cadu
% Marito
% Agus
% Euge
% -------------------------------
% 10-Abril-2019, Universidad Politecnica de Catalunya
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%profile on

% ------------------------
% ****************
% INPUTS
% ****************

% YOUNG's MODULUS
% ---------------
YOUNG_M = 2.00E+11 ;

% Poisson's coefficient
% -----------------------
POISSON = 0.26 ;


% Isotropic modulus
% ---------------------------
K =[0];%YOUNG_M/8;
%K = [ 0.0 YOUNG_M/16 YOUNG_M/8 YOUNG_M/4 YOUNG_M/2 ];

% Kinematic modulus
% ---------------------------
%HMod =[0];
HMod = [ 0.0 YOUNG_M/16 YOUNG_M/8 YOUNG_M/4 YOUNG_M/2 ];

% Modulus for Exponential Hardening
% ---------------------------
DeltaMod = [200.0];
%DeltaMod = [0.0 150 400 800 2000];

% Yield stress
% ------------
YIELD_STRESS = 2.50E+08 ;


% SOFTENING/HARDENING TYPE
% ------------------------
HARDTYPE = 'LINEAR' ; %{PERFECT,LINEAR,EXPONENTIAL}

% VISCOUS/INVISCID
% ------------------------
VISCOUS = 'NO' ;

% Viscous coefficient ----
% ------------------------
eta = [5e10] ;
%eta = [1 2e10 7e10 15e10];

% TimeTotal (initial = 0) ----
% ------------------------
TimeTotal =[10];
%TimeTotal =[100 10 5 2 ];

% Integration coefficient v (for mid-point rule)
% ------------------------
v = 1 ;

% Applied Stress
% ------------------------
nloadstates = 3;
SIGMA = zeros(nloadstates,1);
sigma = 5.0E+08;
SIGMA = [sigma
        -sigma
        sigma];

% Number of time increments for each load state
% --------------------------------------- 
istep=30;

% ------------------------
% ****************
% FUNCTION CALLS
% ****************- 
switch  HARDTYPE
    case 'PERFECT'
        hard_type = 0  ;
    case 'LINEAR'
        hard_type = 1  ;
    case 'EXPONENTIAL'
        hard_type = 2  ;
    otherwise
        hard_type = 0  ;
end

switch  VISCOUS
    case 'NO'
        visc = 0  ;
    case 'YES'
        visc = 1  ;
    otherwise
        visc = 0  ;
end





%%
%%%% SINGLE SIMULATION ------------------------------------------------..s
%matprop=[YOUNG_M,YIELD_STRESS,hard_type,K,HMod, DeltaMod,visc,eta,POISSON];
%Ce=elastic_tensor(matprop);
%[strain_vec,sigma_vec,TIME,dev_sigma_vec]=PlasticityMainJ2(matprop,Ce,STRAIN,TimeTotal,istep);
%%%% SINGLE SIMULATION ------------------------------------------------..e


%%
%%%% MULTI SIMULATION -------------------------------------------------..s
%For Plotting Purposes and comparison of the parameters
%It is needed to change the values in the vectors input of matprop and 
%the function PlasticityMain depending which value want to compare
for L=1:4
    
    matprop=[YOUNG_M,YIELD_STRESS,hard_type,K,HMod(L), DeltaMod(1),visc,eta(1),POISSON];
    Ce=elastic_tensor(matprop);
    STRAIN = jStrain(YOUNG_M,SIGMA,istep,POISSON);
    [strain_vec,sigma_vec,TIME,dev_sigma_vec]=PlasticityMainJ2(matprop,Ce,STRAIN,TimeTotal,istep);
    


    %Colors
    %[0.00 0.45 0.74]blue [0.85 0.33 0.10]orange [0.93 0.69 0.13]yellow [0.64 0.08 0.18]red
    
    
%     % ONE PLOT SIGMA
    figure(1)
    hold on
    if L==1
        plot(strain_vec(1,:),sigma_vec(1,:),'-*','LineWidth',2, 'color',[0 0 0]);
    elseif L==2
        plot(strain_vec(1,:),sigma_vec(1,:),'-*','LineWidth',2, 'color',[0.00 0.45 0.74]);
    elseif L==3
        plot(strain_vec(1,:),sigma_vec(1,:),'-*','LineWidth',2, 'color',[0.85 0.33 0.10]);
    elseif L==4
        plot(strain_vec(1,:),sigma_vec(1,:),'-*','LineWidth',2, 'color',[0.64 0.08 0.18]);
    else
        plot(strain_vec(1,:),sigma_vec(1,:),'-*','LineWidth',2, 'color',[0.93 0.69 0.13]);
    end
        
    xlabel('Strain $\varepsilon_{11}$','Interpreter','latex') 
    ylabel('Stress $\sigma_{11}$','Interpreter','latex') 
    title('Perfect Plas $\nu$ Comp: $\sigma_{11}$ - $\varepsilon_{11}$ Plot','Interpreter','latex') 
    %l=legend('$\eta=1.0$','$\eta=2.0\cdot(10^{10})$','$\eta=7.0\cdot(10^{10})$','$\eta=15\cdot(10^{10})$','Interpreter','latex','Location','northwest');
    %l=legend('$\dot{\varepsilon}=1.3\cdot(10^{-4})$','$\dot{\varepsilon}=1.3\cdot(10^{-3})$','$\dot{\varepsilon}=2.6\cdot(10^{-3})$','$\dot{\varepsilon}=6.5\cdot(10^{-3})$','Interpreter','latex','Location','northwest');
    %l=legend('$K=0$','$K=1.25\cdot(10^{10})$','$K=2.5\cdot(10^{10})$','$K=5\cdot(10^{10})$','$K=1\cdot(10^{11})$','Interpreter','latex','Location','northwest');
    %l=legend('$H=0.0$','$H=1.3\cdot(10^{10})$','$H=2.5\cdot(10^{10})$','$H=5.0\cdot(10^{10})$','$H=1.0\cdot(10^{11})$','Interpreter','latex','Location','northwest');
    %l=legend('$\delta=0$','$\delta=150$','$\delta=500$','$\delta=1500$','$\delta=3000$','Interpreter','latex','Location','northwest');
    set(gca,'FontSize',18);
    %set(gca,'ytick',10^8*([ -3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]))
    grid(gca,'minor')
    grid on

    % ONE PLOT DEVIATORIC SIGMA
    figure(2)
    hold on
    if L==1
        plot(strain_vec(1,:),dev_sigma_vec(1,:),'-*','LineWidth',2, 'color',[0 0 0]);
    elseif L==2
        plot(strain_vec(1,:),dev_sigma_vec(1,:),'-*','LineWidth',2, 'color',[0.00 0.45 0.74]);
    elseif L==3
        plot(strain_vec(1,:),dev_sigma_vec(1,:),'-*','LineWidth',2, 'color',[0.85 0.33 0.10]);
    elseif L==4
        plot(strain_vec(1,:),dev_sigma_vec(1,:),'-*','LineWidth',2, 'color',[0.64 0.08 0.18]);
    else
        plot(strain_vec(1,:),dev_sigma_vec(1,:),'-*','LineWidth',2, 'color',[0.93 0.69 0.13]);
    end
        
    xlabel('Strain $\varepsilon_{11}$','Interpreter','latex') 
    ylabel('Deviatoric Stress $dev(\sigma_{11}$)','Interpreter','latex') 
    title('Perfect Plasticity $\nu$ Comp: $dev(\sigma_{11}$) - $\varepsilon_{11}$ Plot','Interpreter','latex') 
    %l=legend('$\eta=1.0$','$\eta=2.0\cdot(10^{10})$','$\eta=7.0\cdot(10^{10})$','$\eta=15\cdot(10^{10})$','Interpreter','latex','Location','northwest');
    %l=legend('$\dot{\varepsilon}=1.3\cdot(10^{-4})$','$\dot{\varepsilon}=1.3\cdot(10^{-3})$','$\dot{\varepsilon}=2.6\cdot(10^{-3})$','$\dot{\varepsilon}=6.5\cdot(10^{-3})$','Interpreter','latex','Location','northwest');
    %l=legend('$K=0$','$K=1.25\cdot(10^{10})$','$K=2.5\cdot(10^{10})$','$K=5\cdot(10^{10})$','$K=1\cdot(10^{11})$','Interpreter','latex','Location','northwest');
    %l=legend('$H=0.0$','$H=1.3\cdot(10^{10})$','$H=2.5\cdot(10^{10})$','$H=5.0\cdot(10^{10})$','$H=1.0\cdot(10^{11})$','Interpreter','latex','Location','northwest');
    %l=legend('$\delta=0$','$\delta=150$','$\delta=500$','$\delta=1500$','$\delta=3000$','Interpreter','latex','Location','northwest');
    set(gca,'FontSize',18);
    %set(gca,'ytick',10^8*([ -3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]))
    grid(gca,'minor')
    grid on

    
%      figure(3)
%      hold on
%      if L==1
%         plot(TIME,sigma_vec,'-*','LineWidth',2, 'color',[0 0 0]);
%     elseif L==2
%         plot(TIME,sigma_vec,'-*','LineWidth',2, 'color',[0.00 0.45 0.74]);
% 
%     elseif L==3
%         plot(TIME,sigma_vec,'-*','LineWidth',2, 'color',[0.85 0.33 0.10]);
% 
%     elseif L==4
%         plot(TIME,sigma_vec,'-*','LineWidth',2, 'color',[0.64 0.08 0.18]);
%      end
%      title('Time - Stress Comparison','Interpreter','latex');
%      xlabel('Time [sec]','Interpreter','latex') ;
%      ylabel('Stress $\sigma$','Interpreter','latex');
%      l=legend('$\eta=1.0$','$\eta=2.0\cdot(10^{10})$','$\eta=7.0\cdot(10^{10})$','$\eta=15\cdot(10^{10})$','Interpreter','latex','Location','southwest');
%      set(gca,'FontSize',18);
%      grid(gca,'minor')
%      grid on
    
end
%%%% MULTI SIMULATION %%%% -------------------------------------------------..e





