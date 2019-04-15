function strain = iStrain(YOUNG_M,SIGMA,istep)

strain =[0] ;
a=0;
i=1;
j=2;

STRAIN=(SIGMA/YOUNG_M);

for i=1:size(SIGMA)
    if STRAIN(i)>=0
        while round(a,8)<=round(STRAIN(i),8)
            strain_a=STRAIN(i)/istep;
            strain(j)=strain(j-1)+strain_a;
            a=a+strain_a;
            j=j+1;
        end
    else
        while round(a,8)>=round(STRAIN(i),8)
            strain_a=STRAIN(i)/istep;
            strain(j)=strain(j-1)+strain_a;
            a=a+strain_a;
            j=j+1;
        end
    end
end
strain=strain';

   
    