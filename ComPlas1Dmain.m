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



  