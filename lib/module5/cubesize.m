function [edgeind,faceind,cnr]=cubesize(b,cubesz,nx,ny,nz)
%indiziert die Kanten eines kleinen oder großen Würfel
%kleiner Würfel cubesz=0: wie in Task5
%großer Würdel cubesz=1 doppelte Größe wie in Task 5

if b==0 && cubesz==1
%größerer Würfel

cnr=9*7*5; %315
edgeind=zeros(cnr,54);
faceind=zeros(cnr,6);

c2=1;
   for kc=1:nz-2        
    for jc=1:ny-2
        for ic=1:nx-2
            
[edgeindt, faceindt]= edgesofcube (ic, jc, kc, nx, ny, nz,cubesz,b);

 edgeind(c2,:)=edgeindt;
 faceind(c2,:)=faceindt;
 c2=c2+1;
        end 
    end 
   end     

else
 
cnr=480;
 c1=1;
edgeind=zeros(cnr,12);
faceind=zeros(cnr,6);

for kc=1:nz-1       %b=0 -1
    for jc=1:ny-1
        for ic=1:nx-1
            
[edgeindt, faceindt]= edgesofcube (ic, jc, kc, nx, ny, nz,cubesz,b);

 edgeind(c1,:)=edgeindt;
 faceind(c1,:)=faceindt;
 c1=c1+1;
        end 
    end 
end
  
        
end