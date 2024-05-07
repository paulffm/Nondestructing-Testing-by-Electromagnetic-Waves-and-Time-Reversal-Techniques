

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
%Für Vorwärtsrichtung
%eps: b=0 eps=eps0 b=1 eps=eps0*80


    load('UCubesizeposs.mat')


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
                                            
% Knu = curlcurl(mesh,1/mu0*(mesh.elem2regi==LB));   % (alternative)       % [1/H] : curl-reluctance-curl matrix


%möglichen Cubes: cubes in x-Richtung xc=floor(wx/wcube)
                 %cubes in y-Richtung yc=floor(wy/wcube)
                 %cubes in z-Richtung zc=floor(wz/wcube)
                 %Anzahl cubes cnr=cx*yc*zc Hier 480

cnr=9*7*5; %315
edgeind=zeros(cnr,54);
%faceind=zeros(cnr,6);


b=0;
%Alle möglichen Kanten und Flächen-Indizes der Cubes
cc=1;
cubesize=1;

for kc=1:nz-2        
    for jc=1:ny-2
        for ic=1:nx-2
            
[edgeindt, faceindt]= edgesofcube (ic, jc, kc, nx, ny, nz,cubesize,b);

 edgeind(cc,:)=edgeindt;
 %faceind(cc,:)=faceindt;
 cc=cc+1;
        end 
    end 
end

%nicht nicht x=3 nicht 23, sondern 23+2*11*9=221
%edgeind(221,:)

%% E. Boundary conditions (and removing "hair edges")
idxpe = setdiff(find(mesh.ds), ...
  idxboundary(mesh,'left/right/front/back/top','pe'));  
 %Position wurde so gewählt, dass in jeder Raumrichtung ein Würfel zwischen
 %dem ersten Würfel und diesem liegt
 %1+1+1=3 in x-Richtung
 %3+2*nx=3+9+9=21 in y-Richtung
 %21+nx*ny*2=21+9*7*2=147
 
 %pos=147;
 pos=74;
 %pos=randi(315); 
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


ia=3;
  for ja=3:7
        na=ia+(ja-1)*nx;
   if ja<7
       
        XA(na)=1;
        XA(na+1)=1;

        XA(na+np)=1;
        XA(na+1+np)=1;
         XA(na+2+np)=1;
 
   else
       
        XA(na)=1;
        XA(na+1)=1;
       
   end
   
  end


% Kanten Elektrode B
ib=7;
  for jb=3:7
 
nb=ib+(jb-1)*nx;

    if jb<7
        XB(nb)=1;
        XB(nb+1)=1;
 
        XB(nb+np)=1;
        XB(nb+1+np)=1;
        XB(nb+2+np)=1;
        
    else 
        XB(nb)=1;
        XB(nb+1)=1;
    end 
    
  end

XAp=XA(idxpe,:);
XBp=XB(idxpe,:);

 %Max Time Step: CFL
 delta_x=min(hmesh)/2;
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
gauss=(1/sqrt(2*pi)*sigma)*exp(-0.5*(t-t0).^2 ./(sigma^2));

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

   
       
