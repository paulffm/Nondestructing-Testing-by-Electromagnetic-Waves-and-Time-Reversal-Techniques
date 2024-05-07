% Box Resonator
% Forschungspraxis 1 und Forschungspraxis 2+
% TU Darmstadt - Institut fuer Theorie Elektromagnetischer Felder
% Prof. Dr.-Ing. Herbert De Gersem

%% A. Data
plotflag = 1;      % []                                                    : toggling plotting on/off

%Insgesamt Größe in x,y,z Richtung
wx = 0.200; wy = 0.160; wz = 0.120;                                        % [m]   : domain sizes

%Material: Permittivität, Permeabilität
eps0 = 8.85418781762e-12;                                                  % [F/m] : permittivity
mu0 = 4*pi*1e-7;                                                           % [H/m] : permeability

%% B. Analytical calculation of the eigenvalues
% to be computed

%% C. Create model
LB = 1;  % [@]   : label
%Schrittweite in x,y,z Richtung: Entfernung Punkte
hmesh = [ 0.040 ; 0.020 ; 0.050 ];    % [m]   : maximal mesh sizes

%erstellt Brick, Würfel mit Ursprung [ -wx/2 -wy/2 -wz/2 ]
%und Kantenlängen [ wx wy wz ]
shp = Brick([ -wx/2 -wy/2 -wz/2 ],[ wx wy wz ]);                           % []    : brick

%erstellt aus den Bricks das und Seitengröße das Mesh
mesh = HexMesh('box','shapes',{shp},LB,hmesh,plotflag);                    % []    : hexahedral mesh

%Plottet Mesh
if plotflag
  figure(1); clf; plotregion(mesh,[],[],[],[],struct('marker','-'));
end

%% D. Create matrices
%Anzahl der Gitterpunkte in Raumrichtung und Mz=nx*ny
[np,nx,ny,nz,nxy,stride] = sizes(mesh);                                    % [#]   : sizes

%Materialmatrix Meps Meps=DAd*Deps*DS
Meps = edgemass(mesh,eps0*(mesh.elem2regi==LB));                           % [F]   : capacitance matrix

%Materialmatrix Mnu=DSd*Dmnu*DA
Mnu = facemass(mesh,1/mu0*(mesh.elem2regi==LB));                           % [H]   : reluctance matrix

%Curl-Curl Matrix
Knu = mesh.C'*Mnu*mesh.C;                                                  % [1/H] : curl-reluctance-curl matrix


%% E. Boundary conditions (and removing "hair edges")
%Elektrische RB: Tang. Kanten werden gelöscht
idxpe = idxboundary(mesh,'left/right/front/back/bottom/top','pe');         % [@]   : indices of all tangential edges on the front, back, under and upper boundaries

%Geisterkanten werden zusätzlich gelöscht
idxdof = setdiff(find(mesh.ds),idxpe);                                     % [@]   : indices of the degrees of freedom

%Auswertung der Matrizen an diesen übrigen Kanten
Knured = Knu(idxdof,idxdof);                                               % [1/H] : curl-reluctance-curl matrix (with boundary conditions)
Mepsred = Meps(idxdof,idxdof);                                             % [F]   : capacitance matrix (with boundary conditions)

%% F. Eigenvalue solution
%Löst Eigenwertproblem A*V = B*V*D.
%mit A=C*Mnu*C und B=Meps V=Eigenvektorn D=Eigenvalues=Kreisfrequenz^2
[V,D] = eig(full(Knured),full(Mepsred));                                   % [V,(rad/s)^2] : eigenvectors and eigenvalues

%Berechnung Eigenfrequenzen aus Kreisfrequenzen^2
f = sqrt(diag(D))/(2*pi);                                                  % [Hz]  : eigenfrequencies

%Plot Eigenfrequenzen
figure(2); clf; plot(real(f),imag(f),'x'); axis([ 0 max(real(f)) -max(abs(imag(f))) max(abs(imag(f))) ]);
xlabel('eigenfrequencies (Hz)'); ylabel('imaginary parts (Hz)');
