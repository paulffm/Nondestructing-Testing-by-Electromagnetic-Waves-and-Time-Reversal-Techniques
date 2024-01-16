

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
%c Material zu eps0=eps0*80
c=0;

%Größe Elektroden
 %b=0 Standard
 %b=1 beide klein
 %b=2 B klein
 %b=3 A klein
b=0;

if c==0
    
    switch b
        case 0    
            load('Uposs.mat')
        case 1
            load('UsEposs.mat')
        case 2
            load('UsBposs.mat')
         case 3
            load('UsAposs.mat')
    end
else 
    load('Uepsposs.mat')
    eps0=eps0*80;
    b=0;
end

%compressed data
[U,S,V] = svd(uBfposs);
uBfposscomp=U(:,1:2)*S(1:2,1:2)*V(:,1:2)';


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
                                            

%möglichen Cubes: Anzahl cubes: ci Cubes in jede Raumrichtung cnr=cx*yc*zc=480 
                 
%Größe Würfel
%b=0 wie in Task5, b=1 doppelt so groß
cubesize=0;

    
cnr=480;
edgeind=zeros(cnr,12);
faceind=zeros(cnr,6);

%Alle möglichen Kanten und Flächen-Indizes der Cubes
%cc Nummerierung nach kanonischen Indize
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


%% E. Boundary conditions (and removing "hair edges")
idxpe = setdiff(find(mesh.ds), ...
  idxboundary(mesh,'left/right/front/back/top','pe'));  

 pos=183; %nach kanonischer Indizierung: Würfelnummer nach
 %Aufgabenstellung
 %pos=90;
 
 %pos=randi(480); 
for ii=1:size(edgeind,2)
   
idxpe(idxpe==edgeind(pos,ii))=[];
end

%Nicht-triviale Kanten=automatisch magn. RB
idxde=find(diag(Mnu)); 

Deps = eps0*mesh.dad(idxpe)./mesh.ds(idxpe);   

Dnu = 1/mu0*mesh.dsd(idxde)./mesh.da(idxde);  

Cdp = mesh.C(idxde,idxpe);   

%Elektroden A,B in x,y-Fläche 
[XA,XB]=sizeelectrodes(b,nx,ny,nz);

%Auswertung Elektroden an nicht Geisterkanten und mit elektr. RB
XAp=XA(idxpe,:);
XBp=XB(idxpe,:);

 %Max Time Step: CFL
 delta_x=min(hmesh)/2;
 deltaTmaxCFL = sqrt( eps0*mu0 ) * sqrt( 1/( 1/(delta_x^2)+ 1/(delta_x^2) + 1/(delta_x^2) ));
 
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
gauss=(1/sqrt(2*pi)*sigma)*exp(-0.5*(t-t0).^2 ./(sigma^2));

%Gauß moduliert
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
 
   %Position des Cubes finden
   %Überprüfung der Signale mit Korrelation
   %sortieren nach größter Korrelation
   %löschen aller außer der drei größten KorrKoeff
   %Ausgabe der KorrKoeff mit Nummerierung
   
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
   Vb(1:(length(Vb)-5))=[];
   Indiceb(1:(length(Indiceb)-5))=[];
   %[Vb,Indiceb]=max(Rb);
   

   
   
   %Kanten des Cubes
   %edgeindf=edgeind(Indicef,:);
   %edgeindb=edgeind(Indiceb,:);
   
   fprintf('Die Position des Würfels ist: %d \n', pos);
   fprintf('Die im Vorwärtsschritt bestimmte Position ist: %d \n',Indicef);
   fprintf('KorrKoff Vorwärtsschritt: %4.3f \n',Vf);
   fprintf('Die nach der Komprimierung bestimmte Position ist: %d \n',Indicefcomp);
   fprintf('KorrKoff Komprimierung: %4.3f \n',Vfcomp);
   fprintf('Die nach dem Time Reversal bestimmte Position ist: %d \n', Indiceb);
   fprintf('KorrKoff Rückwärtsschritt: %4.3f \n',Vb);
 
   

   

       
