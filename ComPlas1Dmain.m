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

<<<<<<< HEAD
% Young's modulus
YOUNG_M = 2.00E+11;

% Poisson's coefficient
POISSON = 0.26 ;

%Kinematic hardening modulus
H = -0.25;

%Isotropic hardening modulus
K = -0.5;

%Yield limit
YIELD_STRESS = 2.50E+08 ;

%Isotropic hardening type (LINEAR-0,EXPONENTIAL-1)
HARDTYPE = 0 ;

% Viscous/Inviscid
VISCOUS = 0 ;

% Viscous coefficient
eta = 1 ;

% TimeTotal (initial = 0
TimeTotal = 10;

% Integration coefficient ALPHA
ALPHA_COEFF = 1 ;

%Maximum load applied
SIGMA = 350;

%Discretization
istep = 10;

%Elastic stress path
strain_path = iStrain(YOUNG_M,SIGMA,istep);



  
=======
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
incrSigma = SIGMA/istep;
incrStrain = iStrain(YOUNG_M,SIGMA,istep);

>>>>>>> Plasticity1
