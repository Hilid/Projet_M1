clear all
close all
clc
graphics_toolkit('fltk')          %affichage gnuplot

%===============================================================================================================
%Constantes physiques
%--------------------
c = 340;
rho = 1.177;       %a  300°K


% Constantes guide
%-----------------
L = 3*0.3; 				% longueur du guide
d = 0.05; 				% diametre du guide


% Constantes Résonateur
%----------------------
Lcav =0.165;			% longueur de la cavité
Lcol =2*0.02;			% longueur du col
Dcav =0.05;			% diametre de la cavité
Dcol =0.02;			% diametre du col du col

Scav = pi*(Dcav/2)^2;
Scol = pi*(Dcol/2)^2;

RN = Dcol / 2;
RC = Dcav / 2;
RT = d / 2;

%Correction de longueur du col (prise dans [A1] appendice B)
L1 = 0.82 * (1 - 1.35*RN/RC + 0.31*(RN/RC)^3) * RN; 
L2 = 0.82 * (1- 0.235 * RN / RT - 1.32*(RN/RT)^2 + 1.54 * (RN/RT)^3 - 0.86*(RN/RT)^4)*RN;
Lcol = Lcol + L1 + L2;



%Base fréquentielle
%------------------
Fmax=2000;
f = 0:1:Fmax;
N = length(f);
w = 2*pi*f;

%=================================================================================================================

%Calcul des coefficients de la matrice
%--------------------------------------
res = ones(1,N);			%matrice de transfert du résonateur
reseau = ones(2,2,N);		%matrice de transfert du réseau composé des éléments de 'config.txt'
admittance_resonateur = ones(1,N);

for x=1:1:N
	reseau(:,:,x) = eye(2);	%initialisation de la matrice par une matrice diago de 1
end


for x=1:1:N
	w = 2*pi* x / N * Fmax ;
	matrice_resonateur = resonateur(w,Lcav,Lcol,Dcav,Dcol,rho,c);
	cellule_reso = guide(w,L,d,rho,c)* matrice_resonateur* guide(w,L,d,rho,c);
	reseau(:,:,x) =reseau(:,:,x) * cellule_reso; 
	admittance_resonateur(x) = matrice_resonateur(2,1);
end

%====================================================================
% Affichage dans le terminal
%====================================================================
disp('===============================================================');
disp(['Paramètres globaux']);
disp('------------------');
disp(['célérité: ',num2str(c)]);
disp(['masse volumique: ',num2str(rho)]);
disp('');
disp(['Paramètres résonateur']);
disp('---------------------');
disp(['Lcav = ' num2str(Lcav)]);
disp(['Lcol = ' num2str(Lcol)]);
disp(['Dcav = ' num2str(Dcav)]);
disp(['Dcol = ' num2str(Dcol)]);
disp('');
disp(['Paramètres guide']);
disp('---------------------');
disp(['Diametre du guide = ' num2str(d)]);
disp(['Longueur entre chaques resonateurs 2L = ' num2str(2*L)]);
disp(['Frequence de la bande de bragg pour c = ' num2str(c) 'm.s^1,   =>    f = c/(4L) = ' num2str(c/(4*L)) ' Hz']);
disp('');
disp('===============================================================');
disp('');

%====================================================================
% Affichage des courbes
%====================================================================

%Affichage de l'équation de dispersion 2cos(Gamma d) = T11 + T12
%--------------------------------------------------------------
cosGammad = ((reseau(1,1,:) + reseau(2,2,:)) /2);
Gd = acos(cosGammad);


figure(1)
subplot(3,1,1)
plot(f,cosGammad,'-b');
hold on
plot(ones(1,N), '-r');
hold on
plot(-ones(1,N), '-r');
axis([0 2000 -10 10])
legend('cos( Gamma d )','y=1','y=-1');
ylabel('cos( Gamma d)');
title("Representation graphique de l'equation de dispersion")
grid minor on

subplot(3,1,2)
plot(f, real(Gd));
ylabel('Re(Gamma d)');
grid minor on

subplot(3,1,3)
plot(f,imag(Gd));
ylabel('Im(Gamma d)');
xlabel('frequence en Hz');
grid minor on




%Affichage de l'impédance caractéristique du réseau Z_reseau
%-----------------------------------------------------
Z_reseau = sqrt(reseau(1,2,:)./reseau(2,1,:));

figure(2)
subplot(4,1,1);
plot(f,log(abs(Z_reseau)), 'r');
%xlabel('frequence en Hz');
ylabel('abs(Zreseau)');
title('Impedance caracteristique du reseau');
grid on


%Affichage du coefficient de reflexion
%--------------------------------------------------------------
S = pi*(d/2)^2;     %section du guide
Zc = rho*c/S;       % Sans perte


% Zc avec pertes
%---------------
%~ Pr = 0.708;                                        % nombre de Prandtl a la pression atmosphérique           a 300°K
%~ mu = 1.85 * 10^-5;                                 % viscosité de l'air		                                a 300°K
%~ gamma = 1.4;                               	      % heat capacity ratio of air
%~ 
%~ ksi = sqrt(Pr);
%~ delta = sqrt(2 * mu / (rho * w));            	  % viscous boundary layer thickness
%~ s = d / 2 / delta;
%~ beta = (1-i)/sqrt(2);
%~ Zc = rho * c / S * (1 + beta /s * (1 - (gamma - 1)/ksi));

R = (reseau(1,1,:) + reseau(1,2,:)./Zc - reseau(2,1,:).*Zc- reseau(2,2,:))./(reseau(1,1,:) + reseau(2,1,:).*Zc + reseau(1,2,:)./Zc + reseau(2,2,:));

figure(2)
subplot(4,1,2);
plot(f,abs(R));
ylabel('abs(R)');
title('Coefficient de reflexion a l entree du reseau')
grid on


%Affichage du coefficient de transmission
%--------------------------------------------------------------
T = 2./(reseau(1,1,:)  + reseau(2,1,:).*Zc + reseau(1,2,:)./Zc + reseau(2,2,:));

figure(2)
subplot(4,1,3);
plot(f,abs(T));
ylabel('abs(T)');
title('Coefficient de transmission a l entree du reseau');
grid on



%Affichage de l'admittance du résonateur
%-----------------------------------------------
figure(2)
subplot(4,1,4);
plot(f,abs(admittance_resonateur));
ylabel('abs(Yreso)');
xlabel('frequence en Hz');
title('admittance du resonateur de Helmholtz');
grid on

