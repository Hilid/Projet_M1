clear all
%close all
clc
graphics_toolkit('fltk')          %affichage gnuplot

%===============================================================================================================
nb_cellule =5;

%Constantes physiques
%--------------------
c = 340;
rho = 1.177;       %a  300°K


% Constantes guide
%-----------------
L = 0.10; 				% longueur du guide
d = 0.05; 				% diametre du guide


% Constantes Résonateur
%----------------------
Lcav =0.16;			% longueur de la cavité
Lcol =0.02;			% longueur du col
Dcav =0.043;			% diametre de la cavité
Dcol =0.02;				% diametre du col

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
f = 0:0.5:Fmax;
N = length(f);
w = 2*pi*f;


%Prise en compte du chaos (changement aléatoire de Lcav avec un écart-type de sigma)
%-----------------------------------------------------------------------------------
sigma =0.1;  %ecart type en mm
vec_L= L+sigma*randn(1,nb_cellule);

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
	for y=1:1:nb_cellule
		reseau(:,:,x) = reseau(:,:,x)* (guide(w,vec_L(y),d,rho,c)*resonateur(w,Lcav,Lcol,Dcav,Dcol,rho,c));
	end
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
disp(['Longueur entre chaques resonateurs L = ' num2str(L) ' m']);
disp(['Frequence de la bande de bragg pour c = ' num2str(c) ' m.s^1,   =>    f = c/(2L) = ' num2str(c/(2*L)) ' Hz']);
disp('');
disp('===============================================================');
disp('');

%====================================================================
% Affichage des courbes
%====================================================================

%Affichage de l'équation de dispersion 2cos(Gamma d) = T11 + T12
%--------------------------------------------------------------
cosnGammad = (A + D) /2;
Gd = acos(cosnGammad)/nb_cellule;


figure(1)
subplot(1,3,1)
plot(cosnGammad,f,'-b');
hold on
plot(ones(1,N),f, '-r');
hold on
plot(-ones(1,N),f, '-r');
axis([ -10 10 0 2000])
%legend('cos( nGamma d )','y=1','y=-1');
xlabel('cos( Gamma d)');
%title("Representation graphique de l'equation de dispersion")
grid minor on

subplot(1,3,2)
plot(real(Gd),f);
xlabel('Re(Gamma d)');
grid minor on

subplot(1,3,3)
plot(imag(Gd),f);
xlabel('Im(Gamma d)');
ylabel('frequence en Hz');
grid minor on


%Affichage de l'impédance caractéristique du réseau Zc
%-----------------------------------------------------
Zr = sqrt(B./C);

%~ figure(2)
%~ subplot(3,1,1);
%~ plot(f,log(abs(Zr)), 'r');
%~ ylabel('log(abs(Zreseau))');
%~ title('Impedance caracteristique du reseau');
%~ grid on


%Affichage du coefficient de reflexion
%--------------------------------------------------------------
S = d^2/4*pi;
Zc = rho*c/S;

Zc_perte = ones(1,1,N);

for x=1:1:N
		w = 2*pi* x / N * Fmax ;
		[Zc_perte(1,1,x) b] = pertes(d,w,rho,c);
end
Zc = Zc_perte;


R = (A + B./Zc - C.*Zc - D)./(A + C.*Zc + B./Zc + D);

figure(2)
hold on
subplot(2,1,1);
hold on
plot(f,20*log10(abs(R)))%,'--r');
%~ axis([0 Fmax 0 1.5]);

ylabel('20log(|R|)');
title("Coefficient de réflexion à l'entrée du réseau")
grid on


%Affichage du coefficient de transmission
%--------------------------------------------------------------
T = 2./(A + C.*Zc + B./Zc + D);

figure(2)
hold on
subplot(2,1,2);
hold on
plot(f,20*log10(abs(T(1,1,:))))%,'--r');
axis([0 Fmax -100 0]);
ylabel('20log(|T|)');
title("Coefficient de transmission à l'entrée du réseau");
grid on



%Affichage de l'admittance du résonateur
%-----------------------------------------------
%~ figure(2)
%~ subplot(3,1,3);
%~ semilogy(f,abs(admittance_resonateur(1,:)));
%~ ylabel('|Yr|');
%~ xlabel('Fréquence en Hz');
%~ title('Admittance du résonateur de Helmholtz');
%~ grid on

legend('Sans desordre','Avec desordre')

print -dsvg chaos_position_grand.svg

%~ %Absorption A=1-abs(T)^2 -abs(R)^2
%~ A=1-abs(T).^2 -abs(R).^2;
%~ 
%~ figure(3)
%~ plot(f,abs(T),'b');
%~ hold on
%~ plot(f,abs(R),'r');
%~ hold on
%~ plot(f,A,'k');
%~ xlim([0 1500]);
%~ legend('Transmission','Reflexion','Absorption');
%~ hold off
%~ title('chaotique')

