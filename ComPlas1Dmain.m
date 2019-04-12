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
POISSON = 0.3 ;

% Plastic modulus
% ---------------------------
K = -1.5 ;

% Kinematic Hardening/softening modulus
% ---------------------------
H = -1.5 ;

% Yield stress
% ------------
YIELD_STRESS = 2.50E+08 ;


% SOFTENING/HARDENING TYPE
% ------------------------
HARDTYPE = 'EXPONENTIAL' ; %{LINEAR,EXPONENTIAL}
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
sigma = 3.50E+08;
SIGMA = [sigma
        -sigma
        sigma];

% Number of time increments for each load state
% --------------------------------------- 
istep=10;
matprop=[YOUNG_M,YIELD_STRESS];

STRAIN = iStrain(YOUNG_M,SIGMA,istep);

sigma_vec=PlasticityMain(matprop,STRAIN,SIGMA,TimeTotal,istep);





