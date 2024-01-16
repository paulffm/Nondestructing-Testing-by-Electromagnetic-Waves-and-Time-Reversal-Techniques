
% Forschungspraxis 1 und Forschungspraxis 2+
% TU Darmstadt - Institut fuer Theorie Elektromagnetischer Felder
% Prof. Dr.-Ing. Herbert De Gersem

%% A. Data
plotflag = 1;                                                              % []    : toggling plotting on/off
wx = 0.200; wy = 0.160; wz = 0.120;                                        % [m]   : domain sizes
eps0 = 8.85418781762e-12;                                                  % [F/m] : permittivity
mu0 = 4*pi*1e-7;                                                           % [H/m] : permeability

%% B. Analytical calculation of the eigenvalues
% anaEW=load('AnalytischEW');
% anaEW=cell2mat(struct2cell(anaEW));

%% C. Create model
LB = 1;                                                                    % [@]   : label

hmesh = [ 0.030 ; 0.030 ; 0.030 ];                                         % [m]   : maximal mesh sizes
shp = Brick([ -wx/2 -wy/2 -wz/2 ],[ wx wy wz ]);                           % []    : brick
mesh = HexMesh('box','shapes',{shp},LB,hmesh,plotflag);                    % []    : hexahedral mesh
if plotflag
  figure(1); clf; plotregion(mesh,[],[],[],[],struct('marker','-'));
end

%% D. Create matrices
[np,nx,ny,nz,nxy,stride] = sizes(mesh);                                    % [#]   : sizes
Meps = edgemass(mesh,eps0*(mesh.elem2regi==LB));                           % [F]   : capacitance matrix
Mnu = facemass(mesh,1/mu0*(mesh.elem2regi==LB));                           % [H]   : reluctance matrix
Knu = mesh.C'*Mnu*mesh.C;                                                  % [1/H] : curl-reluctance-curl matrix
% Knu = curlcurl(mesh,1/mu0*(mesh.elem2regi==LB));   % (alternative)       % [1/H] : curl-reluctance-curl matrix

%% E. Boundary conditions (and removing "hair edges")
idxpe = idxboundary(mesh,'left/right/front/back/bottom/top','pe');         % [@]   : indices of all tangential edges on the front, back, under and upper boundaries
idxdof = setdiff(find(mesh.ds),idxpe);                                     % [@]   : indices of the degrees of freedom

%Werte an existierenden Kanten aus
%sonst konvergiert Solver später nicht
%[exedg,ghostedg]=getGhostEdges(nx,ny,nz);

%Auswertung an nicht Geisterkanten+elektrische RB
%sonst teilen durch 0
Knured = Knu(idxdof,idxdof);                                                 % [1/H] : curl-reluctance-curl matrix (with boundary conditions)
%Mepsred = Meps(exedg,exedg); % [F]   : capacitance matrix (with boundary conditions)

%Auswertung an nicht Geisterkanten+elektrische RB
%sonst teilen durch 0
Mepsred=Meps(idxdof,idxdof);


%Stromanregung
%b=0 Leiterschleife; b=1 einzelne Anregung an einer Kante; b=2 Anregung an
%einem Punkt in alle Richtung
%b=3 zufällige Anregung an 1/3 aller Kanten
 b=3; 

%Amplitude der Stromstärke
 ampl=1;
%jbow Stromverteilung, indices= an welchen Kanten allokiert
[jbow, indices]=currentdistribution(b,ampl,nx,ny,nz,idxdof);
 
 %E soll nur an existierenden Kanten ausgewertet werden+elektr RB
 ebow=zeros(size(idxdof,1),3000);
 
 %Impedanzvektor
 zimp=zeros(3000,1);
 
 %Frequenzvektor
 freqv=zeros(3000,1);
 
 for i=1:length(freqv)
     
 f=10.^(8+i*0.001); 
 freqv(i)=f;
 
 omega=2*pi*f;

 j=sqrt(-1);
 
 mkap=0;
 
 %Systemmatrix 
 Ared = Knured+j*omega*mkap-omega^2*Mepsred;

 %rechte Seite nur existEdges+RB
 rhsred=-j*omega*jbow(idxdof);
 x=Ared\rhsred;
 ebow(:,i)=x;           %x ist LSG des GLS und damit elektrisches Feld
 
 
 %zimp(i,1)=(x(9)+x(13)+x(45)+x(46))/10;
 %Z=U/I an Kanten, an denen Strom eingespeist wird
 
 %Da jbow für die Berechnung von x nur an exist. Kanten ausgwertet wurde,
 %haben sich die Indizes der Allokation des Stromes verschoben
 %Zur Berechnung der Impedanz an den Kanten der Allokation des Stroms
 %muss also die Spannung an diesen "verschobenen Kanten ausgewertet werden
 %Daher mit find nicht Nulleinträge finden, und x an diesen Kanten
 %auswerten
 
 zimp(i,1)=sum(x(find(jbow(idxdof))))./ampl;                                             
 
 end
 
%Betrag,Winkel Impedanz 
zabs=abs(zimp);
zangle=angle(zimp);

%Plot Betrag Z
figure (1); clf; loglog(freqv,zabs);
xlabel('Frequenz (Hz)'); ylabel('Magnitude of Z');
title('Bodeplot');
grid on;

%Plot Winkel Z +/- 90 Grad
figure (2); clf; semilogx(freqv,zangle);
xlabel('Frequenz (Hz)'); ylabel('Phase of Z');
title('Bodeplot');
 grid on;

 %Berechnung Eigenfrequenzen
 freso=freqv(islocalmax(zabs));
 
 %Vergleich mit analytischen EW zum Test für mich
 fa=load('AnalytischeEW');
 fa=cell2mat(struct2cell(fa));
 fnum=freso(1:15,:);
 fa=fa(1:15,1);
 relerror=abs(fnum-fa)./fa;
 sumrel=sum(relerror);
 
