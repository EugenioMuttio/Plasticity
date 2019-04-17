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

% Plastic modulus
% ---------------------------
K =YOUNG_M/4  ;

% Kinematic Hardening/softening modulus
% ---------------------------
H = 0.0 ;

% Yield stress
% ------------
YIELD_STRESS = 20 ;


% SOFTENING/HARDENING TYPE
% ------------------------
HARDTYPE = 'PERFECT' ; %{PERFECT,LINEAR,EXPONENTIAL}

% VISCOUS/INVISCID
% ------------------------
VISCOUS = 'NO' ;

% Viscous coefficient ----
% ------------------------
eta = 1 ;

% TimeTotal (initial = 0) ----
% ------------------------
TimeTotal = 10 ;

% Integration coefficient v (for mid-point rule)
% ------------------------
v = 1 ;

% Applied Stress
% ------------------------
nloadstates = 3;
SIGMA = zeros(nloadstates,1);
sigma = 35;
SIGMA = [sigma
        -sigma
        sigma];

% Number of time increments for each load state
% --------------------------------------- 
istep=50;

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


matprop=[YOUNG_M,YIELD_STRESS,hard_type,K,H];

STRAIN = iStrain(YOUNG_M,SIGMA,istep);

[strain_vec,sigma_vec]=PlasticityMain(matprop,STRAIN,SIGMA,TimeTotal,istep);

plot(strain_vec,sigma_vec,'-o');

%%TEST
nstrain=size(strain_vec);
strstr=zeros(nstrain(1),2);
strstr(:,1)=strain_vec;
strstr(:,2)=sigma_vec;
grid on;





