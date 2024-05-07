function [jbow, indices]= currentdistribution(i,ampl,nx,ny,nz,idxdof)
np=nx*ny*nz;
jbow=zeros(3*np,1);

%Stromanregung
%i=0 Leiterschleife
%i=1 Anregung einer Kante
%i=2 Anregung alle von einem Punkt abgehenden Kanten
%Zufallige Anregung an 1/3 aller Kanten
%i=3 zufällig, ob +- Polarisation?
%i=4 +

  switch i
      case 0
   n1=1+(3-1)+(3-1)*nx+(2-1)*nx*ny;         %i=3, j=3 k=2 
   n2=1+(3-1)+(4-1)*nx+(2-1)*nx*ny;          %i=3, j=4 k=2
   n3=1+(3-1)+(3-1)*nx+(2-1)*nx*ny+np;      %i=3, j=3 k=2 +np=
   n4=1+(4-1)+(3-1)*nx+(2-1)*nx*ny+np;      %i=4  j=3 k=2 +np=162 

 jbow(n1)=ampl;   
 jbow(n2)=-ampl;  
 jbow(n3)=-ampl; 
 jbow(n4)=ampl;  
 
 indices=[n1 n2 n3 n4];
 
      case 1
        
        n1=3+nx+nx*ny;
        jbow(n1)=ampl;
        indices=n1;
           
      
      case 2
              
              
        %falsche Zuordnung ich erreiche ja alles damit
        n1=1+(3-1)+(3-1)*nx+(2-1)*nx*ny; %i=j=3; k=2
        
        
        jbow(n1)=ampl;
        jbow(n1+np)=ampl;
        jbow(n1+2*np)=ampl;        
        indices=[n1 n1+np n1+2*np];
        
          
      case 3  
    length3=floor(length(idxdof)/3);
   auswahl=randperm(length(idxdof));
    indices=idxdof(auswahl(1:length3));
    for i=1:length(indices)
        a=randsample([-1, 1],1,true);
    jbow(indices(i))=a*ampl;
    
    end
       
      case 4
          length3=floor(length(idxdof)/3);
   auswahl=randperm(length(idxdof));
    indices=idxdof(auswahl(1:length3));

    jbow(indices)=ampl;
          
          
  end   

            
            
   
        
