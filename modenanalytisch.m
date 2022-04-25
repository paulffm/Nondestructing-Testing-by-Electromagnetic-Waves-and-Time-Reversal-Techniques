fate=zeros(96,1);
fatm=zeros(80,1);
mnpte=strings(96,1);
mnptm=strings(80,1);

eps0 = 8.85418781762e-12;                                                  % [F/m] : permittivity
mu0 = 4*pi*1e-7;                                                           % [H/m] : permeability
wx = 0.200; wy = 0.160; wz = 0.120;            
%TE p>=1
ite=1;
itm=1;
for mte=0:4
    for nte=0:4
        for pte=1:4
 
            if mte==0 && nte==0
          
               
            else
                
fate(ite)=(1/sqrt(eps0*mu0))*sqrt((mte/(2*wx)).^2 + (nte/(2*wy)).^2 + (pte/(2*wz)).^2);

      mnpte(ite,1)=['TE' num2str(mte) num2str(nte) num2str(pte)];
       ite=ite+1;

            end
  
        end
    end 
end

%TM m>=1 n>=1;

for mtm=1:4
    for ntm=1:4
        for ptm=0:4
           
fatm(itm)=(1/sqrt(eps0*mu0))*sqrt((mtm/(2*wx)).^2 + (ntm/(2*wy)).^2 + (ptm/(2*wz)).^2);            
       
 mnptm(itm,1)=['TM' num2str(mtm) num2str(ntm) num2str(ptm)];

itm=itm+1;
        end 
    end 
end

[fate,Ite]=sort(fate);
[fatm,Itm]=sort(fatm);

mnpte=mnpte(Ite,:);
mnptm=mnptm(Itm,:);

fa=[fate; fatm];
mnp=[mnpte; mnptm];

[fa,It]=sort(fa);
mnp=mnp(It,:);

cut=123; 

fa=fa(1:cut);
mnp=mnp(1:cut,:);


    
