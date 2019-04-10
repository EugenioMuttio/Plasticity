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
% Hardening/softening modulus
% ---------------------------
HARDSOFT_MOD = -1.5 ;
% Yield stress
% ------------
YIELD_STRESS = 2.50E+08 ;
% Problem type  TP = {'PLANE STRESS','PLANE STRAIN','3D'}
% ------------------------ = 1            =2         =3
% ------------
ntype= 2 ;
% Model    PTC = {'SYMMETRIC','TENSION','NON-SYMMETRIC'} ;
%                     = 1         = 2         = 3
% ---------------------------------------------------
MDtype =3;
% Ratio compression strength / tension strength
% ---------------------------------------------
n = 4 ;
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
TimeTotal = 10 ; ;
% Integration coefficient ALPHA
% ------------------------
ALPHA_COEFF = 1 ;
% Points ---------------------------
% ----------------------------------
nloadstates = 3 ;
SIGMAP = zeros(nloadstates,2) ;
SIGMAP(1,:) =[3.75E+08 0];
SIGMAP(2,:) =[-3.00E+08 0];
SIGMAP(3,:) =[5.75E+08 0];
% Number of time increments for each load state
% --------------------------------------- 
istep = 10*ones(1,nloadstates) ;
