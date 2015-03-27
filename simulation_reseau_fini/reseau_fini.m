clear all
close all
clc
graphics_toolkit('gnuplot')          %affichage gnuplot

%===============================================================================================================
nb_cellule =20;

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
Fmax=750;
f = 0:0.5:Fmax;
N = length(f);
w = 2*pi*f;

%=================================================================================================================

%Calcul des coefficients de la matrice
%--------------------------------------
reseau = ones(2,2,N);		%matrice de transfert du réseau composé des éléments de 'config.txt'
admittance_resonateur = ones(1,N);


for x=1:1:N
	matrice_resonateur = resonateur(w(x),Lcav,Lcol,Dcav,Dcol,rho,c);
	cellule = guide(w(x),L,d,rho,c) * matrice_resonateur;

	reseau(:,:,x) = eye(2);
	for y=1:1:nb_cellule
		reseau(:,:,x) = reseau(:,:,x)*cellule;
	end
	
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
legend('cos( nGamma d )','y=1','y=-1');
xlabel('cos( Gamma d)');
title("Representation graphique de l'equation de dispersion")
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

figure(2)
subplot(4,1,1);
plot(f,log(abs(Zr)), 'r');
ylabel('log(abs(Zreseau))');
title('Impedance caracteristique du reseau');
grid on


%Affichage du coefficient de reflexion
%--------------------------------------------------------------
S = d^2/4*pi;
Zc = rho*c/S;

Zc_perte = ones(1,1,N);

for x=1:1:N
		[Zc_perte(1,1,x) b] = pertes(d,w(x),rho,c);
end
Zc = Zc_perte;


R = (A + B./Zc - C.*Zc - D)./(A + C.*Zc + B./Zc + D);

figure(2)
subplot(4,1,2);
plot(f,abs(R));
axis([0 Fmax 0 1.5]);

ylabel('abs(R)');
title('Coefficient de reflexion a l entree du reseau')
grid on


%Affichage du coefficient de transmission
%--------------------------------------------------------------
T = 2./(A + C.*Zc + B./Zc + D);

figure(2)
subplot(4,1,3);
plot(f,abs(T(1,1,:)));
axis([0 2000 0 1.5]);
ylabel('abs(T)');
title('Coefficient de transmission a l entree du reseau');
grid on



%Affichage de l'admittance du résonateur
%-----------------------------------------------
figure(2)
subplot(4,1,4);
semilogy(f,abs(admittance_resonateur(1,:)));
ylabel('abs(Yreso)');
xlabel('frequence en Hz');
title('admittance du resonateur de Helmholtz');
grid on

%Absorption A=1-abs(T)^2 -abs(R)^2
A=1-abs(T).^2 -abs(R).^2;

figure(3)
plot(f,abs(T),'b');
hold on
plot(f,abs(R),'r');
hold on
plot(f,A,'k');
xlim([0 Fmax]);
legend('Transmission','Reflexion','Absorption');
hold off
title('periodique')
