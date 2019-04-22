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
YOUNG_M = 20000 ;

% Poisson's coefficient
% -----------------------
POISSON = 0.3 ;

% Isotropic modulus
% ---------------------------
K =0.0%YOUNG_M/4;

% Kinematic modulus
% ---------------------------
HMod =0;%YOUNG_M/4;

% Modulus for Exponential Hardening
% ---------------------------
DeltaMod = 5000.0;

% Yield stress
% ------------
YIELD_STRESS = 20 ;


% SOFTENING/HARDENING TYPE
% ------------------------
HARDTYPE = 'PERFECT' ; %{PERFECT,LINEAR,EXPONENTIAL}

% VISCOUS/INVISCID
% ------------------------
VISCOUS = 'YES' ;

% Viscous coefficient ----
% ------------------------
eta = 5e3 ;

% TimeTotal (initial = 0) ----
% ------------------------
TimeTotal = 100 ;

% Integration coefficient v (for mid-point rule)
% ------------------------
v = 1 ;

% Applied Stress
% ------------------------
nloadstates = 3;
SIGMA = zeros(nloadstates,1);
sigma = 50;
SIGMA = [sigma
        -sigma
        sigma];

% Number of time increments for each load state
% --------------------------------------- 
istep=40;

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

matprop=[YOUNG_M,YIELD_STRESS,hard_type,K,HMod, DeltaMod,visc,eta];

STRAIN = iStrain(YOUNG_M,SIGMA,istep);

[strain_vec,sigma_vec,TIME]=PlasticityMain(matprop,STRAIN,SIGMA,TimeTotal,istep);

figure(1)
hold on
plot(strain_vec,sigma_vec,'-o');

% figure(2)
% hold on
% plot(TIME,sigma_vec,'-o');


%%TEST
nstrain=size(strain_vec);
strstr=zeros(nstrain(1),2);
strstr(:,1)=strain_vec;
strstr(:,2)=sigma_vec;
grid on;





