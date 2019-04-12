function sigma_vec = PlasticityMain(matprop,STRAIN,SIGMA,TimeTotal,istep)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returning Map Algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%5 (paths)
delta_t=TimeTotal/istep/5;

% Strain Rate
eps_rate=[];
sigma_vec=zeros(size(STRAIN));
for i=1:size(STRAIN)-1
    i=i+1;
    eps_n=STRAIN(i-1);
    eps_n1=STRAIN(i);
    eps_rate(i)=(STRAIN(i)-STRAIN(i-1))/delta_t;
    
    sigma_vec(i)=maps_plas(matprop,eps_rate,eps_n,eps_n1);
end





    