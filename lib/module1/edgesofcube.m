function [edgeind, faceind]= edgesofcube (i, j, k, nx, ny, nz,cubesz,b)
%größer machen wäre
%dann bspw 24 Kanten
if b==0 && cubesz==1
 
 N0=1+ (i-1) +nx* (j-1) +ny*nx* (k-1) ;
 Np=nx*ny*nz;
edgeind=zeros(54,1);

%x
edgeind (1) =N0 ; %bottom (x)
edgeind(2)=N0+nx; %(x)

edgeind (3) =N0+1 ; %bottom (x)
edgeind(4)=N0+1+nx; %(x)

edgeind (5) =N0+2*nx ; %bottom (x)
edgeind(6)=N0+2*nx+1; %(x)

%y
edgeind (7) =N0 +Np;
edgeind(8)=N0+1+Np; 
edgeind (9) =N0+2;

edgeind(10)=N0+nx+Np; 
edgeind(11)=N0+nx+1+Np; 
edgeind(12)=N0+nx+2+Np; 



edgeind (13:24) =edgeind (1:12) +ny*nx;
edgeind(25:36)=edgeind(13:24)+ny*nx;

edgeind(37)=N0+2*Np; %z
edgeind(38)=N0+1+2*Np;
edgeind(39)=N0+2+2*Np; %z

edgeind(40)=N0+nx+2*Np;
edgeind(41)=N0+nx+1+2*Np;
edgeind(42)=N0+nx+2+2*Np;

edgeind(43)=N0+2*nx+2*Np;
edgeind(44)=N0+2*nx+1+2*Np;
edgeind(45)=N0+2*nx+2+2*Np;

edgeind(46:54)=edgeind(37:45)+ny*nx;

faceind=zeros(6,1);
faceind(1) =N0 ; %x
faceind(2)=N0+Np; %y
faceind(3)=N0+2*Np;
faceind(4)=N0+nx*ny+2*Np;%z
faceind(5)=N0+ny+Np;%y
faceind(6)=N0+1;%x       
        

else
    
    
N0=1+ (i-1) +nx* (j-1) +ny*nx* (k-1) ;
Np=nx*ny*nz;
edgeind=zeros(12,1);
edgeind (1) =N0 ; %bottom (x)
edgeind(2)=N0+nx; %(x)
edgeind (3) =N0 +Np;%(y)
edgeind(4)=N0+1+Np; %(x)

edgeind (5:8) =edgeind (1:4) +ny*nx; %top

edgeind(9)=N0+2*Np; %z
edgeind(10)=N0+1+2*Np;
edgeind(11)=N0+nx+2*Np;
edgeind(12)=N0+nx+1+2*Np;

faceind=zeros(6,1);
faceind (1) =N0 ; %x
faceind (2)=N0+Np; %y
faceind(3)=N0+2*Np;
faceind(4)=N0+nx*ny+2*Np;%z
faceind(5)=N0+ny+Np;%y
faceind(6)=N0+1;%x
          
        
end
    
