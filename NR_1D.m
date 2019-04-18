function gamma_n1 = NR_1D(func,dfdxi,ftrial,mat_p)

%Material Properties
xi_n=mat_p(1);
E=mat_p(2);
H=mat_p(3);
K=mat_p(4);
delta=mat_p(5);
sigma_y=mat_p(6);
sigma_inf=mat_p(7);

% Variables Initialization
k=0;
gamma_n1=0;
g_res=1;

%Tolerance
tol=1*10^(-12);
%Max Iterations
Max_Iter=50;

while (abs(g_res)>tol && k<Max_Iter)
    xi_inc=xi_n+gamma_n1;
    % Residual G 
    g_res=ftrial-gamma_n1*(E+H)-(func(xi_inc,sigma_inf,sigma_y,delta,K)-...
        func(xi_n,sigma_inf,sigma_y,delta,K));
    %Derivative of the residual
    dgdxi=-(E+dfdxi(xi_inc,delta,K)+H);
    
    delta_gamma=-g_res/dgdxi;

    gamma_n1=delta_gamma+gamma_n1;
    k=k+1;
    %fprintf("Iteration %d - gamma= %e - g(xi)=%e\n",k,gamma_n1,g_res);
end
fprintf("FOUND: Iteration %d - gamma= %e- g(xi)=%e \n",k,gamma_n1,g_res);



