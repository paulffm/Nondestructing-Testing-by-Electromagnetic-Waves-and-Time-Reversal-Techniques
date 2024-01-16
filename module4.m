

% Box Resonator
% Forschungspraxis 1 und Forschungspraxis 2+
% TU Darmstadt - Institut fuer Theorie Elektromagnetischer Felder
% Prof. Dr.-Ing. Herbert De Gersem


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
if plotflag
  figure(1); clf; plotregion(mesh,[],[],[],[],struct('marker','-'));
end

%% D. Create matrices
[np,nx,ny,nz,nxy,stride] = sizes(mesh);                                    % [#]   : sizes
Meps = edgemass(mesh,eps0*(mesh.elem2regi==LB));                           % [F]   : capacitance matrix
Mnu = facemass(mesh,1/mu0*(mesh.elem2regi==LB));                           % [H]   : reluctance matrix
                                            
% Knu = curlcurl(mesh,1/mu0*(mesh.elem2regi==LB));   % (alternative)       % [1/H] : curl-reluctance-curl matrix

%% E. Boundary conditions (and removing "hair edges")
idxpe = setdiff(find(mesh.ds), ...
  idxboundary(mesh,'left/right/front/back/top','pe'));  

%Task4 Magnetische RB am Boden, sonst elektrische
%was ist Boden? x,y Fläche
%Wie implementiere ich elektrische? mit idxpe automatisch
%->idxpe alle tangentialen Kanten am Rand

%Nicht-triviale Kanten=automatisch magn. RB, wieso eig.?
idxde=find(diag(Mnu)); 

%Meps=DAd*Deps*DS Deps=flächengemitt Permittivität 
Deps = eps0*mesh.dad(idxpe)./mesh.ds(idxpe);   

%Mnu-1=DSd*Dmnu-1*DA Dmnu=längengemitt Permeabilität
Dnu = 1/mu0*mesh.dsd(idxde)./mesh.da(idxde);  

%Curl Matrix
Cdp = mesh.C(idxde,idxpe);   

 %Elektroden A,B in x,y-Fläche:
 
 %Größe der Elektroden
 %b=0 Standard
 %b=1 beide klein
 %b=2 B klein
 %b=3 A klein
  b=0;
  
 [XA,XB]=sizeelectrodes(b,nx,ny,nz);


XAp=XA(idxpe,:);
XBp=XB(idxpe,:);
 
 
 %Max Time Step: CFL
 delta_x=0.02; 
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

 for nf=1:length(t) 
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

 for nb=1:length(t)
  hb = hb-dt*Dnu.*(Cdp*eb);
  eb = eb+dt*(Cdp'*hb-XBp*iB(nb))./Deps;
  uAb(nb,1) = XAp'*eb;
  uBb(nb,1)=XBp'*eb;
 end
  
 
 %Plot Forward
  figure(2);
  tiledlayout(2,2)
 
  nexttile
  plot(t./1e-9,uAf*1e3); xlabel('time (ns)'); ylabel('voltage (mV)');
  title("voltage electrode A forward")
  
  nexttile
  plot(t./1e-9,uBf*1e3); xlabel('time (ns)'); ylabel('voltage (mV)');
  title("voltage electrode B forward")
  
  nexttile([1 2])
  plot(t./1e-9,iF*1e6); xlabel('time (ns)'); ylabel('current (uA)');
  title("current forward")



  %Plot backwards
  figure(3)
  tiledlayout(2,2)
 
  nexttile
  plot(t./1e-9,uAb*1e3); xlabel('time (ns)'); ylabel('voltage (mV)');
  title("voltage electrode A backward")
  
  nexttile
  plot(t./1e-9,uBb*1e3); xlabel('time (ns)'); ylabel('voltage (mV)');
  title("voltage electrode B backward")
  
  nexttile([1 2])
  plot(t./1e-9,iB*1e3); xlabel('time (ns)'); ylabel('current (mA)');
  title("current backward")
 
 %Lege Signal aus Rückwärts über Signal Vorwärts an Elektrode A
 uAbflip=flip(uAb);
 figure(4)
 plot(t./1e-9,uAf); xlabel('time (ns)'); ylabel('voltage (V)');
 xlabel('time (ns)'); ylabel('voltage (V)');
 title('Vergleich Spannung an A vorwärts und rückwärts');
 hold on
 %standarf 1.4e-5
 plot(t./1e-9,uAbflip*1.4e-5); xlabel('time (ns)'); ylabel('voltage (V)');
 legend({'U_A vorwärts','U_A rückwärts'},'Location','EastOutside')

 
