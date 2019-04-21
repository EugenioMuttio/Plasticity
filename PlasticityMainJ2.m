function [strain_vec,sigma_vec,TIME]= PlasticityMainJ2(matprop,Ce,STRAIN,TimeTotal,istep)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returning Map Algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%5 (paths)
delta_t=TimeTotal/istep/5;

% Strain Rate
eps_rate=zeros(size(STRAIN));

%viscous
visc=matprop(7);

eps_pvec=zeros(size(STRAIN));
xi_vec=zeros(size(STRAIN));%scalar
xibar_vec=zeros(size(STRAIN));%vector
q_vec=zeros(size(STRAIN));%scalar
qbar_vec=zeros(size(STRAIN));%vector
gamma_vec=zeros(size(STRAIN));

sigma_vec=zeros(size(STRAIN));
strain_vec=zeros(size(STRAIN));
TIME=zeros(1,size(STRAIN,2)); % initialize time variable

for i=1:size(STRAIN,2)-1
    i=i+1;
    TIME(i)=TIME(i-1)+delta_t;
    %strains
    eps_n=STRAIN(:,i-1);
    eps_n1=STRAIN(:,i);
    eps_rate(:,i)=(STRAIN(:,i)-STRAIN(:,i-1))/delta_t;
    %internal variable: plastic strains
    eps_p_n=eps_pvec(:,i-1);
    eps_p_n1=eps_pvec(:,i);
    %internal variable: xi
    xi_n=xi_vec(:,i-1);
    xi_n1=xi_vec(:,i);
    %internal variable: xibar
    xibar_n=xibar_vec(:,i-1);
    xibar_n1=xibar_vec(:,i);
    %hardening variable: q
    q_n=q_vec(:,i-1);
    q_n1=q_vec(:,i);
    %hardening variable: qbar
    qbar_n=qbar_vec(:,i-1);
    qbar_n1=qbar_vec(:,i);
    
    %internal variable: gamma
    gamma_n=gamma_vec(:,i-1);
    gamma_n1=gamma_vec(:,i);
    
    int_vars_nn1=[eps_n eps_n1 eps_p_n eps_p_n1 xi_n xi_n1 xibar_n xibar_n1 gamma_n gamma_n1 q_n q_n1 qbar_n qbar_n1];
    
    if visc==0
        [sigma_vec(:,i),int_vars_nn1]=maps_plasJ2(matprop,Ce,sigma_vec(:,i-1),eps_rate(:,i),int_vars_nn1,delta_t);
    else
        [sigma_vec(:,i),int_vars_nn1]=maps_visplas(matprop,Ce,sigma_vec(:,i-1),eps_rate(:,i),int_vars_nn1,delta_t);
    end
    
    %strains
    eps_pvec(:,i)=int_vars_nn1(:,4);
    xi_vec(:,i)=int_vars_nn1(:,6);
    xibar_vec(:,i)=int_vars_nn1(:,8);
    gamma_vec(:,i)=int_vars_nn1(:,10);
    q_vec(:,i)=int_vars_nn1(:,12);
    qbar_vec(:,i)=int_vars_nn1(:,14);
    
    strain_vec(:,i)=STRAIN(:,i);
end






    