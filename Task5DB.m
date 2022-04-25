

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
mu0 = 4*pi*1e-7;                                                           % [H/m] : permeabilit



%% C. Create model
LB = 1;                                                                    % [@]   : label
hmesh = [ 0.040 ; 0.040 ; 0.040 ]/2;                                       % [m]   : maximal mesh sizes
shp = Brick([ -wx/2 -wy/2 -wz/2 ],[ wx wy wz ]);                           % []    : brick
mesh = HexMesh('box','shapes',{shp},LB,hmesh,plotflag);                    % []    : hexahedral mesh
%if plotflag
 % figure(1); clf; plotregion(mesh,[],[],[],[],struct('marker','-'));
%end

%% D. Create matrices
[np,nx,ny,nz,nxy,stride] = sizes(mesh);                                    % [#]   : sizes
Meps = edgemass(mesh,eps0*(mesh.elem2regi==LB));                           % [F]   : capacitance matrix
Mnu = facemass(mesh,1/mu0*(mesh.elem2regi==LB));                           % [H]   : reluctance matrix
                                            
% Knu = curlcurl(mesh,1/mu0*(mesh.elem2regi==LB));   % (alternative)       % [1/H] : curl-reluctance-curl matrix


%möglichen Cubes: poss=(nx-1)*(ny-1)*(nz-1)=480
cubesz=1;
b=0;

cnr=480;
edgeind=zeros(cnr,12);
faceind=zeros(cnr,6);

cnr=9*7*5;
edgeind=zeros(cnr,54);
faceind=zeros(cnr,6);


cc=1;
for kc=1:nz-2%b=0 -1 überall %b=1 -2 überall
    for jc=1:ny-2
        for ic=1:nx-2
            
[edgeindt, faceindt]= edgesofcube (ic, jc, kc, nx, ny, nz,cubesz,b);
 edgeind(cc,:)=edgeindt;
 faceind(cc,:)=faceindt;
 cc=cc+1;
        end 
    end 
end



 %Elektroden A,B in x,y-Fläche 
 %Größe Elektroden
 %b=0 Standard
 %b=1 beide klein
 %b=2 B klein
 %b=3 A klein

 [XA,XB]=sizeelectrodes(b,nx,ny,nz);
  
 %Max Time Step: CFL
 delta_x=min(hmesh)/2;
 deltaTmaxCFL = sqrt( eps0*mu0 ) * sqrt( 1/( 1/(delta_x^2)+ 1/(delta_x^2) + 1/(delta_x^2) )) ;

 
%Time-Step
dt=deltaTmaxCFL*0.9; 
tend=20e-9; 
t = [0:dt:tend]';

%zu 1 geändert anstatt 1e-
%Signal
ihat = 1e-6;                %Amplitude des Signals
f=5e9;                      %Frequenz
y1 = ihat*sin(2*pi*f*t);

%0.5 geändert anstatt 5
%verschobener Gaußpuls
t0=1e-9;   
sigma=0.5e-9; 
gauss=(1/sqrt(2*pi)*sigma)*exp(-0.5*(t-t0).^2 ./(sigma^2));

%Signal zur Anregung vorwärts
iF=gauss.*y1;

%RB von Gitter
idxpe = setdiff(find(mesh.ds), ...
  idxboundary(mesh,'left/right/front/back/top','pe'));  

%Nicht-triviale Kanten=automatisch magn. RB
idxde=find(diag(Mnu));   

Deps = eps0*mesh.dad(idxpe)./mesh.ds(idxpe);   

Dnu = 1/mu0*mesh.dsd(idxde)./mesh.da(idxde);  

Cdp = mesh.C(idxde,idxpe); 

uAfposs= zeros(length(t),cnr);     %Spannung Elektrode A
uBfposs= zeros(length(t),cnr);     %Spannung Elektrode B

uAbposs= zeros(length(t),cnr);     %Spannung Elektrode A
uBbposs= zeros(length(t),cnr);     %Spannung Elektrode B

for iii=1:cnr

%Nötig, damit immer nur Rand und aktuelle Position des Würfels RB haben
idxpei=idxpe;
    
Depsi =Deps;  

Cdpi =Cdp; 

%RB von Cube: PEC
for ii=1:size(edgeind,2)
    
idxpei(idxpei==edgeind(iii,ii))=[];

end 

Depsi = eps0*mesh.dad(idxpei)./mesh.ds(idxpei);    

Cdpi = mesh.C(idxde,idxpei);  

XAp=XA(idxpei,:);
XBp=XB(idxpei,:);



%Leapfrog forward
ef= zeros(length(idxpei),1);  %E-Feld                                              
hf= zeros(length(idxde),1);  %H-Feld


 for nf=1:length(t)-1 
  hf = hf-dt*Dnu.*(Cdpi*ef);
  ef = ef+dt*(Cdpi'*hf-XAp*iF(nf))./Depsi;
  uAfposs(nf,iii) = XAp'*ef;
  uBfposs(nf,iii)=XBp'*ef;
 end
 
 %Leapfrog backward
eb= zeros(length(idxpei),1);  %E-Feld                                             
hb= zeros(length(idxde),1);  %H-Feld

%Signal backwards
%resistance=uAf./iF;
%iB=uBf./resistance;
iB=flip(uBfposs);

 %eb=ef;
 %hb=hf;
 %uBb=

 for nb=1:length(t)-1 
  hb = hb-dt*Dnu.*(Cdpi*eb);
  eb = eb+dt*(Cdpi'*hb-XBp*iB(nb))./Depsi;
  uAbposs(nb,iii) = XAp'*eb;
  uBbposs(nb,iii)=XBp'*eb;
  end
 
end


switch b
        case 0    
            save('UCubesizeposs','uAfposs','uAbposs','uBfposs','uBbposs')
        case 1
            save('UsEposs','uAfposs','uAbposs','uBfposs','uBbposs')
        case 2
           save('UsBposs','uAfposs','uAbposs','uBfposs','uBbposs')
         case 3
            save('UsAposs','uAfposs','uAbposs','uBfposs','uBbposs')
end







 
 


