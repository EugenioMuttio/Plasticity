function [strain_vec,sigma_vec]= PlasticityMain(matprop,STRAIN,SIGMA,TimeTotal,istep)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returning Map Algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%5 (paths)
delta_t=TimeTotal/istep/5;

% Strain Rate
eps_rate=[];

%viscous
visc=matprop(7);

eps_pvec=zeros(size(STRAIN));
xi_vec=zeros(size(STRAIN));
xibar_vec=zeros(size(STRAIN));
gamma_vec=zeros(size(STRAIN));

sigma_vec=zeros(size(STRAIN));
strain_vec=zeros(size(STRAIN));

for i=1:size(STRAIN)-1
    i=i+1;
    %strains
    eps_n=STRAIN(i-1);
    eps_n1=STRAIN(i);
    eps_rate(i)=(STRAIN(i)-STRAIN(i-1))/delta_t;
    %internal variable: plastic strains
    eps_p_n=eps_pvec(i-1);
    eps_p_n1=eps_pvec(i);
    %internal variable: xi
    xi_n=xi_vec(i-1);
    xi_n1=xi_vec(i);
    %internal variable: xibar
    xibar_n=xibar_vec(i-1);
    xibar_n1=xibar_vec(i);
    
    %internal variable: gamma
    gamma_n=gamma_vec(i-1);
    gamma_n1=gamma_vec(i);
    
    int_vars_nn1=[eps_n eps_n1 eps_p_n eps_p_n1 xi_n xi_n1 xibar_n xibar_n1 gamma_n gamma_n1];
    
    if visc==0
        [sigma_vec(i),int_vars_nn1]=maps_plas(matprop,sigma_vec(i-1),eps_rate(i),int_vars_nn1,delta_t);
    else
        [sigma_vec(i),int_vars_nn1]=maps_visplas(matprop,sigma_vec(i-1),eps_rate(i),int_vars_nn1,delta_t);
    end
    
    %strains
    eps_pvec(i)=int_vars_nn1(4);
    xi_vec(i)=int_vars_nn1(6);
    xibar_vec(i)=int_vars_nn1(8);
    gamma_vec(i)=int_vars_nn1(10);
    
    strain_vec(i)=STRAIN(i);

    
end





    