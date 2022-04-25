

% Box Resonator
% Forschungspraxis 1 und Forschungspraxis 2+
% TU Darmstadt - Institut fuer Theorie Elektromagnetischer Felder
% Prof. Dr.-Ing. Herbert De Gersem


%% A. Data
plotflag = 1;                                                              % []    : toggling plotting on/off
wx = 0.200; wy = 0.160; wz = 0.120;                                        % [m]   : domain sizes
eps0 = 8.85418781762e-12;                                                  % [F/m] : permittivity
mu0 = 4*pi*1e-7;     % [H/m] : permeabilit



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
%Kanten, die nicht tangential am Rand sind, und keine Geisterkanten
%sonst teilt man durch 0
idxpe = setdiff(find(mesh.ds), ...
  idxboundary(mesh,'left/right/front/back/bottom/top','pe'));  

%Ausgwertet an Nicht-Nulleinträgen
%sonst teilt man durch 0
idxde=find(diag(Mnu)); 

%Meps=DAd*Deps*DS Deps=flächengemitt Permittivität (normalerweise)
Deps = eps0*mesh.dad(idxpe)./mesh.ds(idxpe);   

%Mnu-1=DSd*Dmnu-1*DA Dmnu=längengemitt Permeabilität (normalerweise)
Dnu = 1/mu0*mesh.dsd(idxde)./mesh.da(idxde);  

%Curl Matrix
Cdp = mesh.C(idxde,idxpe);   

%Stromanregung
%b=0 Leiterschleife; b=1 einzelne Anregung an einer Kante; b=2 Anregung an
%einem Punkt in alle Richtung
%b=3 zufällige Anregung an 1/3 aller Kanten
 b=3;
 ampl=1;
 
 [X, indices]= currentdistribution(b,ampl,nx,ny,nz,idxpe);
 
 %Stromverteilung ausgwertet auf nicht tang. Kanten
 Xp=X(idxpe,:);
 
 %Max Time Step:
 delta_x=min(hmesh); %aus Mesh /2
 deltaTmaxCFL = sqrt( eps0*mu0 ) * sqrt( 1/( 1/delta_x^2 + 1/delta_x^2 + 1/delta_x^2 ) );

 
%Time-Step
dt=deltaTmaxCFL; 
tend=100e-9; %100e-9
t = [0:dt:tend]';

%Frequenz-Step
N=size(t,1);
Fs=1/dt; %sampling rate
dF = Fs/N; %frequency steps                   
F = (0:floor(N/2)-1)*dF;


%Chirp
ihat = 1e-6;  %i in A
t1= 20e-9;  %ab hier keine weitere Frequenzerhöhung
f1 = 1e9;   %Signal soll bei t=0 Frequenz von f1 besitzen
f2=5e9;     %Signal soll bei t1 Frequenz von f2 besitzen
           %ab dann keine weitere Erhöhung der Frequenz 
           
%Definition Chirp als Sinus-Signal          
y1 = ihat*chirp(t,f1,t1,f2,'linear',-90);

%verschobener Gaußpuls
t0=10e-9;                    %Verschiebung
sigma=5e-9;                  %Standardabweichung
gauss= (1/(sqrt(pi*2)*sigma))*exp(-0.5*(t-t0).^2 ./(sigma^2));

%Signal zur Anregung
i=gauss.*y1;

e = zeros(length(idxpe),1);  %E-Feld                                              % [V]   : grid electric voltages
h= zeros(length(idxde),1);   %H-Feld
u= zeros(length(t),1);       %Spannung U

%Leapfrog Algo nach Skript
 for n=1:length(t)-1 
  h = h-dt*Dnu.*(Cdp*e);
  e = e+dt*(Cdp'*h-Xp*i(n))./Deps;
  u(n,1) = Xp'*e;
 end

 
%Plot Anregung und Spannung im Zeitbereich
if plotflag
  figure(2); clf;
  %/1e-9 *1e6
  subplot(211); plot(t,i); xlabel('time (s)'); ylabel('current (A)');
  title("current")
  subplot(212); plot(t,u); xlabel('time (s)'); ylabel('voltage (V)');
  title("voltage")

end

%Berechnung mit FFT aus Signal im ZB, Signal im FB
%für Spannung und Strom

uf=fft(u);
numel_fft = numel(u); 
y2 = uf(1:floor(N/2));      % -> only positive part is considered
y2=abs(y2);

If=fft(i);
i2=If(1:floor(N/2));
i2=abs(i2);


z2=y2./ihat;

%Plot u im Frequenzbereich
figure(5);
tiledlayout(2,1)
nexttile
plot(F,i2*1e6);xlabel('frequency (Hz)'); ylabel('current (uA)')
xlim([1e9, 5e9])
title("Absolutbetrag des Stroms im FB")


nexttile
plot(F,y2*1e3);xlabel('frequency (Hz)'); ylabel('voltage (mV)')
xlim([1e9, 5e9])
title("Absolutbetrag der Spannung im FB")


figure(6);
plot(F,z2);xlabel('frequency (Hz)'); ylabel('Z(f) (Ohm)')
xlim([1e9, 5e9])
title("Absolutbetrag der Impedanz im FB")



%Berechnung Resonanzfrequenzen
freso=F(islocalmax(y2));

%Überprüfung, ob Resonanzfrequenzen bei Impedanz mit Freq übereinstimmen (nur für mich)
freso2=F(islocalmax(z2));

%Vergleich Analytisch zu numerisch
fa=load('AnalytischeEW');
fa=cell2mat(struct2cell(fa));
 fnum=freso(:,1:15)';
 fa=fa(1:15,1);
 relerror=abs(fnum-fa)./fa;
 sumrel=sum(relerror);
 


%zz=u./i;
%zf=fft(zz);
%zf=If(1:floor(N/2));
%bzf=abs(zf);
%azf=angle(zf);
%plot(F,bzf)
