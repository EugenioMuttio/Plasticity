<<<<<<< HEAD
function [strain_path] = iStrain(YOUNG_M,SIGMA,istep)

strain = zeros(5*istep,1);

delta = SIGMA/(YOUNG_M*istep);

max_strain = SIGMA/YOUNG_M;

for i = 1:size(strain)-1
   if abs(strain(i))<max_strain
       strain(i+1)=strain(i)+delta;
   else
       
   end
    

end

=======
function strain = iStrain(YOUNG_M,SIGMA,istep)

strain = zeros(5*istep+1,1) ;
a=0;
i=1;
j=2;
%k=1;
STRAIN=(SIGMA/YOUNG_M);

for i=1:size(SIGMA)
%     if STRAIN(i)>0
%         k=1;
%     else
%         k=-1;
%     end
    while round(a,8)~=round(STRAIN(i),8)
        strain_a=STRAIN(i)/istep;
        strain(j)=strain(j-1)+strain_a;
        a=a+strain_a;
        j=j+1;
    end
end

    
%     if (SIGMA(i)-a)>0
%         strain_a=(SIGMA(1)/YOUNG_M)*ones(istep+1,1)
%         while 
%         a=SIGMA(i);
%     else
%         strain_a=-(SIGMA(1)/YOUNG_M)*ones(istep+1,1)

    
>>>>>>> Plasticity1
