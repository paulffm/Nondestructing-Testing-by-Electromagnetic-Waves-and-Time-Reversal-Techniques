% Box Resonator
% Forschungspraxis 1 und Forschungspraxis 2+
% TU Darmstadt - Institut fuer Theorie Elektromagnetischer Felder
% Prof. Dr.-Ing. Herbert De Gersem

%% A. Data
plotflag = 1;                                                              % []    : toggling plotting on/off
wx = 0.200; wy = 0.160; wz = 0.120;                                        % [m]   : domain sizes
eps0 = 8.85418781762e-12;                                                  % [F/m] : permittivity
mu0 = 4*pi*1e-7;                                                           % [H/m] : permeability

%% B. Analytical calculation of the eigenvalues
% to be computed     
%initialisiere vektor für Frequenzen TE, TM Moden
fate=zeros(96,1);                                                          %Freq für TE; TM
fatm=zeros(80,1);  

%Init für Speicherung, welche Frequenz, welcher Mode entspricht (mnp Werte)
mnpte=strings(96,1);
mnptm=strings(80,1);                                                       %welche Mode: m,n,p; was 0

%TE,H Ez=0 p>!=1 sonst Hz=0
ite=1;
itm=1;

for mte=0:4
    for nte=0:4
        for pte=1:4 %p>0
 
            if mte==0 && nte==0 %immer nur maximal eins 0 von m,n =0 sonst alles 0 bis auf Hz
          
               
            else
                
fate(ite)=(1/sqrt(eps0*mu0))*sqrt((mte/(2*wx)).^2 + (nte/(2*wy)).^2 + (pte/(2*wz)).^2);
        
      mnpte(ite,1)=['TE' num2str(mte) num2str(nte) num2str(pte)];
       ite=ite+1;

            end
  
        end
    end 
end

%TM,E Hz=0 m>=1 n>=1; sonst alle Komp bis auf Ez =0

for mtm=1:4
    for ntm=1:4
        for ptm=0:4
           
fatm(itm)=(1/sqrt(eps0*mu0))*sqrt((mtm/(2*wx)).^2 + (ntm/(2*wy)).^2 + (ptm/(2*wz)).^2);            
       
 mnptm(itm,1)=['TM' num2str(mtm) num2str(ntm) num2str(ptm)]; %TM 110..

itm=itm+1;
        end 
    end 
end

%sortiere nach Frequenz aufsteigend                     %optionale Sort.
[fate,Ite]=sort(fate);
[fatm,Itm]=sort(fatm);

%sortiere Strings genauso
mnpte=mnpte(Ite,:);
mnptm=mnptm(Itm,:);

%zusammenfügen
fa=[fate; fatm];
mnp=[mnpte; mnptm];

%nur diese Sortierung nötig
%erneut sortieren                                       %wichtige Sort.
[fa,It]=sort(fa);
mnp=mnp(It,:);

cut=15; %15 number of eigenvalues

%schneide ab
%4. und 5. EW sind gleich TE111, TM111: 1,73 Ghz
fa=fa(1:cut);
mnp=mnp(1:cut,:);

%% C. Create model
LB = 1;     % [@]   : label

% 3 Gitter versch. Größen
%i=1 muss kleinstes Gitter sein 

for i=1:3
hmesh = [ 0.020 ; 0.040 ; 0.050 ];                                         % [m]   : maximal mesh sizes
shp = Brick([ -wx/2 -wy/2 -wz/2 ],[ wx wy wz ]);                           % []    : brick
mesh = HexMesh('box','shapes',{shp},LB,hmesh(i)*ones(1,3),plotflag);                    % []    : hexahedral mesh
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
Knured = Knu(idxdof,idxdof);                                               % [1/H] : curl-reluctance-curl matrix (with boundary conditions)
Mepsred = Meps(idxdof,idxdof);                                             % [F]   : capacitance matrix (with boundary conditions)

%% F. Eigenvalue solution
[V,D] = eig(full(Knured),full(Mepsred));                                    % [V,(rad/s)^2] : eigenvectors and eigenvalues

%mit find lösch ich statische Moden und Moden mit Im>0
%Zahlen stehen für verschiedene Gitter
%unsauber mit elseif...

if i==1
fn1 = sqrt(diag(D))/(2*pi); 
fnr1=sort(fn1(find(real(fn1)>10^9)));

elseif i==2
        fn2=sqrt(diag(D))/(2*pi); 
        fnr2=sort(fn2(find(real(fn2)>10^9)));

elseif i==3
            fn3=sqrt(diag(D))/(2*pi); 
            fnr3=sort(fn3(find(real(fn3)>10^9)));

end



end

%größte Größe: wie viele Eigenwerte gibt es
%schreibe alle Frequenzen in diesen Vektor
%Länge ist arbiträr, aber mind. 15 (wenn ich vorher abschneiden würde)
fnr=zeros(size(fnr1,1),3);

%sortierte frequenzen für verschiedene Längen in Spalten schreiben
%fnr entspricht nur diesen dynamischen reellen Eigenwerten
%in den Spalten stehen die EW für verschiedene Gittergrößen
fnr(1:size(fnr1,1),1)=fnr1;
fnr(1:size(fnr2,1),2)=fnr2;
fnr(1:size(fnr3,1),3)=fnr3;

%Frequenzen abschneiden: nur die ersten 15 EW
fnr=fnr(1:cut,:);


%% G. Konvergenz

%fa auf 3 Spalten bringen, damit für alle Gitter Fehler berechnet werden
%kann
fa=fa(1:cut)*ones(1,3);
relerror=abs(fnr-fa)./fa;

%Anpassung des hmesh-Vektors: 15 Zeilen, 3 Spalten
%in allen Zeilen einer Spalte steht Größe des Gitters
%sodass geplottet werden kann; dim(hmeshang)=dim(relerror)
%hmeshang reine Variable zum plotten

hmeshang=zeros(cut,3);
for z=1:cut
    hmeshang(z,:)=hmesh;
end
%Plot: für jeden EW wird rel. Error in Abhängigkeit der Mesh Size geplottet
%also jede "Gerade" entspricht einem der 15 EW und besteht aus 3 Punkten:
% P1=(rel.Error EW(i),mesh(1)) P2=(rel.Error EW(i),mesh(2)) 
% P3=(rel.Error EW(i),mesh(3))
%Nächste Gerade dann entsprechend i+1->nächster EW wird betrachtet
%Knick, da mesh sizes nicht linear (in log) genommen, 0.02->0.04->0.05

%figure (4); clf; loglog(hmeshang',relerror');
%xlabel('mesh size'); ylabel('relative error');
%title('relativer Fehler für die ersten 15 EW');
%legend({'1,199 Ghz','1,456 Ghz','1,561 Ghz','1,732'},'Location','EastOutside')

hmeshang=hmeshang(1:4,:);
relerror=relerror(1:4,:);
figure (4); clf; loglog(hmeshang',relerror');
xlabel('mesh size'); ylabel('relative error');
title('relativer Fehler für die ersten 4 EW');
legend({'1,199 Ghz','1,456 Ghz','1,561 Ghz','1,732 Ghz'},'Location','EastOutside')

figure(5); clf; 
plot(real(fnr(1:3,1)),imag(fnr(1:3,1)),'rx'); 
xlabel('eigenfrequencies (Hz)'); ylabel('imaginary parts (Hz)');
title('Vergleich analytische EW mit den numerischen EW für verschiedene Gittergrößen');
hold on
plot(real(fnr(1:3,2)),imag(fnr(1:3,2)),'b+');
hold on
plot(real(fnr(1:3,3)),imag(fnr(1:3,3)),'c*'); 
hold on
plot(real(fa(1:3,1)),imag(fa(1:3,1)),'ko');
legend({'Mesh size 0.02','Mesh size 0.04','Mesh size 0.05','Analytischen EW'},'Location','EastOutside')

%Konvergenzordnung
ordn=floor((log(relerror(1,3))-log(relerror(1,1)))/(log(hmeshang(1,3))-log(hmeshang(1,1))));
%Fehlerordnung: log(y2)-log(y1)/log(x2)-log(x1) -> ca. 2
fprintf('Konvergenzordnung:%d.\n', ordn)
