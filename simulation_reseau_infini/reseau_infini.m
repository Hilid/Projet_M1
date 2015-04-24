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
L = 0.2; 				% longueur du guide
d = 0.03; 				% diametre du guide


% Constantes Résonateur
%----------------------
Lcav =0.165;			% longueur de la cavité
Lcol =0.02;			% longueur du col
Dcav =0.0425;			% diametre de la cavité
Dcol =0.02;			% diametre du col

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
reseau = ones(2,2,N);		%matrice de transfert du réseau composé des éléments de 'config.txt'
admittance_resonateur = ones(1,N);

for x=1:1:N
	reseau(:,:,x) = eye(2);	%initialisation de la matrice par une matrice identité
end


for x=1:1:N
	w = 2*pi* x / N * Fmax ;
	matrice_resonateur = resonateur(w,Lcav,Lcol,Dcav,Dcol,rho,c);
	reseau(:,:,x) = guide(w,L,d,rho,c) * matrice_resonateur;
	admittance_resonateur(x) = matrice_resonateur(2,1);
end

A = reseau(1,1,:);
B = reseau(1,2,:);
C = reseau(2,1,:);
D = reseau(2,2,:);

%====================================================================
% Affichage dans le terminal
%====================================================================
disp('===============================================================');
disp(['Paramètres globaux']);
disp('------------------');
disp(['célérité: ',num2str(c) ' m.s^-1']);
disp(['masse volumique: ',num2str(rho) ' kg.m^-3']);
disp('');
disp(['Paramètres résonateur']);
disp('---------------------');
disp(['Lcav = ' num2str(Lcav) ' m']);
disp(['Lcol = ' num2str(Lcol) ' m']);
disp(['Dcav = ' num2str(Dcav) ' m']);
disp(['Dcol = ' num2str(Dcol) ' m']);
disp('');
disp(['Paramètres guide']);
disp('---------------------');
disp(['Diametre du guide = ' num2str(d) ' m']);
disp(['Longueur entre chaques resonateurs L = ' num2str(2*L) ' m']);
disp(['Frequence de la bande de bragg pour c = ' num2str(c) ' m.s^1,   =>    f = c/(2L) = ' num2str(c/(2*L)) ' Hz']);
disp('');
disp('===============================================================');
disp('');

%====================================================================
% Affichage des courbes
%====================================================================

%Affichage de l'équation de dispersion 2cos(Gamma d) = T11 + T12
%--------------------------------------------------------------
cosGammad = (A + D) /2;
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


%Affichage de l'impédance caractéristique du réseau Zc
%-----------------------------------------------------
Zc = sqrt(B./C);

figure(2)
subplot(2,1,1)
plot(f,log(abs(Zc)), 'r');
ylabel('log(abs(Zreseau))');
title('Impedance caracteristique du reseau');
grid on


%Affichage de l'admittance du résonateur
%-----------------------------------------------
subplot(2,1,2)
semilogy(f,abs(admittance_resonateur(1,:)));
ylabel('abs(Yreso)');
xlabel('frequence en Hz');
title('admittance du resonateur de Helmholtz');
grid off

