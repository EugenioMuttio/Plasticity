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
K =0;%YOUNG_M/8;


% Kinematic modulus
% ---------------------------
Hmod =[0];
%Hmod = [ 0.0 YOUNG_M/16 YOUNG_M/8 YOUNG_M/4 YOUNG_M/2 ];

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
VISCOUS = 'YES' ;

% Viscous coefficient ----
% ------------------------
%eta = [5e10] ;
eta = [1 2e10 7e10 15e10];

% TimeTotal (initial = 0) ----
% ------------------------
TimeTotal =[10];
TimeTotal =[100 10 5 2 ];

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



STRAIN = iStrain(YOUNG_M,SIGMA,istep);
%%
%%%% SINGLE SIMULATION ------------------------------------------------..s
%matprop=[YOUNG_M,YIELD_STRESS,hard_type,K,Hmod(1), DeltaMod(1),visc,eta(1)];
%[strain_vec,sigma_vec,TIME]=PlasticityMain(matprop,STRAIN,SIGMA,TimeTotal(1),istep);
%%%% SINGLE SIMULATION ------------------------------------------------..e


%%
%%%% MULTI SIMULATION -------------------------------------------------..s
%For Plotting Purposes and comparison of the parameters
%It is needed to change the values in the vectors input of matprop and 
%the function PlasticityMain depending which value want to compare
for L=1:4
    
    matprop=[YOUNG_M,YIELD_STRESS,hard_type,K,Hmod(1), DeltaMod(1),visc,eta(L)];
    [strain_vec,sigma_vec,TIME]=PlasticityMain(matprop,STRAIN,SIGMA,TimeTotal(1),istep);
    

    %Divide the plots to color them different
    nstrain_plot=(size(strain_vec)-1)/5;
    strain_m=zeros(nstrain_plot(1),5);
    sigma_m=zeros(nstrain_plot(1),5);
    K=1;
    for I=1:5
        J=1;
        while J<=nstrain_plot(1)
            strain_m(J,I)=strain_vec(K);
            sigma_m(J,I)=sigma_vec(K);
            J=J+1;
            K=K+1;
        end
        strain_m(J,I)=strain_vec(K);
        sigma_m(J,I)=sigma_vec(K);
    end


    %Colors
    %[0.00 0.45 0.74]blue [0.85 0.33 0.10]orange [0.93 0.69 0.13]yellow [0.64 0.08 0.18]red
    
    
    % ONE PLOT 
    figure(1)
    hold on
    if L==1
        plot(strain_vec(:,1),sigma_vec(:,1),'-*','LineWidth',2, 'color',[0 0 0]);
    elseif L==2
        plot(strain_vec(:,1),sigma_vec(:,1),'-*','LineWidth',2, 'color',[0.00 0.45 0.74]);
    elseif L==3
        plot(strain_vec(:,1),sigma_vec(:,1),'-*','LineWidth',2, 'color',[0.85 0.33 0.10]);
    elseif L==4
        plot(strain_vec(:,1),sigma_vec(:,1),'-*','LineWidth',2, 'color',[0.64 0.08 0.18]);
    else
        plot(strain_vec(:,1),sigma_vec(:,1),'-*','LineWidth',2, 'color',[0.93 0.69 0.13]);
    end
        
    xlabel('Strain $\varepsilon$','Interpreter','latex') 
    ylabel('Stress $\sigma$','Interpreter','latex') 
    title('Nonlinear Iso+Kin $\dot{\varepsilon}$ Comp: $\sigma$ - $\varepsilon$ Plot','Interpreter','latex') 
    %l=legend('$\eta=1.0$','$\eta=2.0\cdot(10^{10})$','$\eta=7.0\cdot(10^{10})$','$\eta=15\cdot(10^{10})$','Interpreter','latex','Location','northwest');
    %l=legend('$\dot{\varepsilon}=1.3\cdot(10^{-4})$','$\dot{\varepsilon}=1.3\cdot(10^{-3})$','$\dot{\varepsilon}=2.6\cdot(10^{-3})$','$\dot{\varepsilon}=6.5\cdot(10^{-3})$','Interpreter','latex','Location','northwest');
    %l=legend('$K=0$','$K=1.25\cdot(10^{10})$','$K=2.5\cdot(10^{10})$','$K=5\cdot(10^{10})$','$K=1\cdot(10^{11})$','Interpreter','latex','Location','northwest');
    l=legend('$H=0.0$','$H=1.3\cdot(10^{10})$','$H=2.5\cdot(10^{10})$','$H=5.0\cdot(10^{10})$','$H=1.0\cdot(10^{11})$','Interpreter','latex','Location','northwest');
    %l=legend('$\delta=0$','$\delta=150$','$\delta=500$','$\delta=1500$','$\delta=3000$','Interpreter','latex','Location','northwest');
    set(gca,'FontSize',18);
    %set(gca,'ytick',10^8*([ -3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]))
    grid(gca,'minor')
    grid on

    %PLOT BY EACH PATH
%     figure(1)
%     hold on
%     plot(strain_m(:,1),sigma_m(:,1),'-*','LineWidth',2, 'color',[0.93 0.69 0.13]);
%     plot(strain_m(:,2),sigma_m(:,2),'-*','LineWidth',2, 'color',[0.85 0.33 0.10]);
%     plot(strain_m(:,3),sigma_m(:,3),'-*','LineWidth',2, 'color',[0.85 0.33 0.10]);
%     plot(strain_m(:,4),sigma_m(:,4),'-*','LineWidth',2, 'color',[0.64 0.08 0.18]);
%     plot(strain_m(:,5),sigma_m(:,5),'-*','LineWidth',2, 'color',[0.64 0.08 0.18]);
%     xlabel('Strain $\varepsilon$','Interpreter','latex') 
%     ylabel('Stress $\sigma$','Interpreter','latex') 
%     title('Nonlinear RD Iso+Kin Plas: $\sigma$ - $\varepsilon$ Plot','Interpreter','latex') 
%     %l=legend('1st: Tensile Load','2nd: Unload/Compression Load',' ','3rd: Unload/Tensile Load','Interpreter','latex','Location','northwest');
%     set(gca,'FontSize',18);
%     grid(gca,'minor')
%     grid on

%      figure(2)
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



