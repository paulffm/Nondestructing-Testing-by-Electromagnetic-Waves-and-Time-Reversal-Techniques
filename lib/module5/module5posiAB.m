

% Box Resonator
% Forschungspraxis 1 und Forschungspraxis 2+
% TU Darmstadt - Institut fuer Theorie Elektromagnetischer Felder
% Prof. Dr.-Ing. Herbert De Gersem

%Script for task5. An obstacle, to be more precise a cube, is introduced within the box resonator.
%A database is build on a signal applied on electrode A and measured at electrode B for varying 
%position of the obstacle. Then signal processing and time reversal is used to determine the 
%position of the box.


%% A. Data
plotflag = 1;                                                              % []    : toggling plotting on/off
wx = 0.200; wy = 0.160; wz = 0.120;                                        % [m]   : domain sizes
eps0 = 8.85418781762e-12;                                                  % [F/m] : permittivity
mu0 = 4*pi*1e-7;   % [H/m] : permeabilit

%Lade Signale an A,B für alle möglichen Positionen von Cube

%ELektroden Angehoben
%d=0 magn. RB am Boden
%d=1 elektr. RB am Boden

d=1;
switch d
    case 0
idxpe = setdiff(find(mesh.ds), ...
  idxboundary(mesh,'left/right/front/back/top','pe'));  
 load('UposiABposs.mat')

    case 1
    idxpe = setdiff(find(mesh.ds), ...
    idxboundary(mesh,'left/right/front/back/bottom/top','pe')); 
    load('UposiABRBposs.mat')
end
   
 

%compressed data
[U,S,V] = svd(uBfposs);
uBfposscomp=U(:,1:20)*S(1:20,1:20)*V(:,1:20)';


%% C. Create model
LB = 1;                                                                    % [@]   : label
hmesh = [ 0.040 ; 0.040 ; 0.040 ]/2;                                       % [m]   : maximal mesh sizes
shp = Brick([ -wx/2 -wy/2 -wz/2 ],[ wx wy wz ]);                           % []    : brick
mesh = HexMesh('box','shapes',{shp},LB,hmesh,plotflag);                    % []    : hexahedral mesh
if plotflag
  figure(1); clf; plotregion(mesh,[],[],[],[],struct('marker','-'));
end


%% D. Create matrices
[np,nx,ny,nz,nxy,stride] = sizes(mesh);                                    % [#]   : sizes
Meps = edgemass(mesh,eps0*(mesh.elem2regi==LB));                           % [F]   : capacitance matrix
Mnu = facemass(mesh,1/mu0*(mesh.elem2regi==LB));                           % [H]   : reluctance matrix
                                            
% Knu = curlcurl(mesh,1/mu0*(mesh.elem2regi==LB));   % (alternative)       % [1/H] : curl-reluctance-curl matrix


%möglichen Cubes: cubes in x-Richtung xc=floor(wx/wcube)
                 %cubes in y-Richtung yc=floor(wy/wcube)
                 %cubes in z-Richtung zc=floor(wz/wcube)
                 %Anzahl cubes cnr=cx*yc*zc Hier 480
%Größe Würfel
%b=0 wie in Task5, b=1 doppelt so groß

    
cnr=480;
edgeind=zeros(cnr,12);
faceind=zeros(cnr,6);

cubesize=0;
%Alle möglichen Kanten und Flächen-Indizes der Cubes
cc=1;

for kc=1:nz-1       %b=0 -1
    for jc=1:ny-1
        for ic=1:nx-1
            
[edgeindt, faceindt]= edgesofcube (ic, jc, kc, nx, ny, nz,cubesize,b);

 edgeind(cc,:)=edgeindt;
 faceind(cc,:)=faceindt;
 cc=cc+1;
        end 
    end 
end

 pos=183;
 %pos=200;
 %pos=randi(480); 
for ii=1:size(edgeind,2)
   
idxpe(idxpe==edgeind(pos,ii))=[];
end

%Nicht-triviale Kanten=automatisch magn. RB
idxde=find(diag(Mnu)); 

Deps = eps0*mesh.dad(idxpe)./mesh.ds(idxpe);   

Dnu = 1/mu0*mesh.dsd(idxde)./mesh.da(idxde);  

Cdp = mesh.C(idxde,idxpe);   

 %Elektroden A,B
 %in x,y-Fläche 
 
 XA=zeros(3*np,1);
 XB=zeros(3*np,1);


% Kanten Elektrode A,B

    %Kanten Elektrode A angehoben
%Elektrode A verschoben

    nxny=3*nx*ny;

    ia=3;
    for ja=3:7
        
        na=ia+(ja-1)*nx;
        if ja<7
       
        XA(na+nxny)=1;
        XA(na+1+nxny)=1;

        XA(na+np+nxny)=1;
        XA(na+1+np+nxny)=1;
         XA(na+2+np+nxny)=1;
        else
     
       
        XA(na+nxny)=1;
        XA(na+1+nxny)=1;
       
        end
   
    end

  
% Kanten Elektrode B angehoben
  ib=7;
  for jb=3:7
 
    nb=ib+(jb-1)*nx;

    if jb<7
        XB(nb+nxny)=1;
        XB(nb+1+nxny)=1;
 
        XB(nb+np+nxny)=1;
        XB(nb+1+np+nxny)=1;
        XB(nb+2+np+nxny)=1;
        
    else 
        XB(nb+nxny)=1;
        XB(nb+1+nxny)=1;
    end 
    
  end

 


XAp=XA(idxpe,:);
XBp=XB(idxpe,:);

 %Max Time Step: CFL
 delta_x=max(hmesh)/2; 
 deltaTmaxCFL = sqrt( eps0*mu0 ) * sqrt( 1/( 1/delta_x^2 + 1/delta_x^2 + 1/delta_x^2 ) );

 
%Time-Step
dt=deltaTmaxCFL*0.9; 
tend=20e-9; 
t = [0:dt:tend]';

%Signal
ihat = 1e-6;                %Amplitude des Signals
f=5e9;                      %Frequenz
y1 = ihat*sin(2*pi*f*t);

%verschobener Gaußpuls
t0=1e-9;   
sigma=0.5e-9; 
gauss= exp(-0.5*(t-t0).^2 ./(sigma^2));

%Signal zur Anregung vorwärts
iF=gauss.*y1;

%Leapfrog forward
ef= zeros(length(idxpe),1);  %E-Feld                                              
hf= zeros(length(idxde),1);  %H-Feld
uAf= zeros(length(t),1);     %Spannung Elektrode A
uBf= zeros(length(t),1);     %Spannung Elektrode B

 for nf=1:length(t)-1 %-1
  hf = hf-dt*Dnu.*(Cdp*ef);
  ef = ef+dt*(Cdp'*hf-XAp*iF(nf))./Deps;
  uAf(nf,1) = XAp'*ef;
  uBf(nf,1)=XBp'*ef;
 end

%Leapfrog backward
eb= zeros(length(idxpe),1);  %E-Feld                                             
hb= zeros(length(idxde),1);  %H-Feld
uAb= zeros(length(t),1);     %Spannung Elektrode A
uBb= zeros(length(t),1);     %Spannung Elektrode B

%Signal backwards
%resistance=uAf./iF;
%iB=uBf./resistance;
iB=flip(uBf);

 %eb=ef;
 %hb=hf;
 %uBb=

 for nb=1:length(t)-1 %-1
  hb = hb-dt*Dnu.*(Cdp*eb);
  eb = eb+dt*(Cdp'*hb-XBp*iB(nb))./Deps;
  uAb(nb,1) = XAp'*eb;
  uBb(nb,1)=XBp'*eb;
 end
 
   %welcher cube
   Rf=corr(uBfposs,uBf);
   [Vf,Indicef]=sort(Rf);
   Vf(1:(length(Vf)-3))=[];
    Indicef(1:(length(Indicef)-3))=[];
    
   %[Vf,Indicef]=max(Rf);
   
   Rfcomp=corr(uBfposscomp,uBf);
   [Vfcomp,Indicefcomp]=sort(Rfcomp);
    Vfcomp(1:(length(Vfcomp)-3))=[];
      Indicefcomp(1:(length(Indicefcomp)-3))=[];
   %[Vfcomp,Indicefcomp]=max(Rfcomp);
   
   
   Rb=corr(uAbposs,uAb);
   [Vb,Indiceb]=sort(Rb);
   Vb(1:(length(Vb)-3))=[];
   Indiceb(1:(length(Indiceb)-3))=[];
   %[Vb,Indiceb]=max(Rb);
   

   
   
   %Kanten des Cubes
   %edgeindf=edgeind(Indicef,:);
  % edgeindb=edgeind(Indiceb,:);
   
   fprintf('Die Position des Würfels ist: %d \n', pos);
   fprintf('Die im Vorwärtsschritt bestimmte Position ist: %d \n',Indicef);
   fprintf('KorrKoff Vorwärtsschritt: %4.3f \n',Vf);
   fprintf('Die nach der Komprimierung bestimmte Position ist: %d \n',Indicefcomp);
   fprintf('KorrKoff Komprimierung: %4.3f \n',Vfcomp);
   fprintf('Die nach dem Time Reversal bestimmte Position ist: %d \n', Indiceb);
   fprintf('KorrKoff Rückwärtsschritt: %4.3f \n',Vb);
 
   

   

       
