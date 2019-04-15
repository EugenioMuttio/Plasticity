function sigma_vec = PlasticityMain(matprop,STRAIN,SIGMA,TimeTotal,istep)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returning Map Algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%5 (paths)
delta_t=TimeTotal/istep/5;

% Strain Rate
eps_rate=[];
eps_p=zeros(size(STRAIN));
sigma_vec=zeros(size(STRAIN));
for i=1:size(STRAIN)-1
    i=i+1;
    eps_n=STRAIN(i-1);
    eps_n1=STRAIN(i);
    eps_p_n=eps_p(i-1);
    eps_p_n1=eps_p(i);
    eps_rate(i)=(STRAIN(i)-STRAIN(i-1))/delta_t;
    
    [eps_p(i),sigma_vec(i)]=maps_plas(matprop,eps_rate(i),eps_n,eps_n1,eps_p_n,eps_p_n1,delta_t);
end





    