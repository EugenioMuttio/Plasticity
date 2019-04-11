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

