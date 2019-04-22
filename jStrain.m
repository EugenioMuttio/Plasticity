function STRAIN = jStrain(YOUNG_M,SIGMA,istep,nu)

strain =[0] ;
a=0;
i=1;
j=2;

LOCSTRAIN=(SIGMA/YOUNG_M);

for i=1:size(SIGMA)
    if LOCSTRAIN(i)>=0
        while round(a,8)<=round(LOCSTRAIN(i),8)
            strain_a=LOCSTRAIN(i)/istep;
            strain(j)=strain(j-1)+strain_a;
            a=a+strain_a;
            j=j+1;
        end
    else
        while round(a,8)>=round(LOCSTRAIN(i),8)
            strain_a=LOCSTRAIN(i)/istep;
            strain(j)=strain(j-1)+strain_a;
            a=a+strain_a;
            j=j+1;
        end
    end
end

STRAIN = [strain; -nu*strain ; -nu*strain ; zeros(1,size(strain,2)); zeros(1,size(strain,2));zeros(1,size(strain,2))];

   
    