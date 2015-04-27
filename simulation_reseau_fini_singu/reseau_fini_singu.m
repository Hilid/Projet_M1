clear all
close all
clc
graphics_toolkit('gnuplot')          %affichage gnuplot

%===============================================================================================================
nb_cellule =9;

%Constantes physiques
c = 340;
rho = 1.177;       %a  300°K

% Constantes guide
L = 0.10; 				% longueur du guide
d = 0.05; 				% diametre du guide

% Constantes Singularité
NumSing = 6;

Lcavs =0.08;			% seul truc qui change
Lcols =0.02;			
Dcavs =0.043;			
Dcols =0.02;				
Scavs = pi*(Dcavs/2)^2;
Scols = pi*(Dcols/2)^2;
RNs = Dcols / 2;
RCs = Dcavs / 2;
RTs = d / 2;
L1s = 0.82 * (1 - 1.35*RNs/RCs + 0.31*(RNs/RCs)^3) * RNs; 
L2s = 0.82 * (1- 0.235 * RNs / RTs - 1.32*(RNs/RTs)^2 + 1.54 * (RNs/RTs)^3 - 0.86*(RNs/RTs)^4)*RNs;
Lcols = Lcols + L1s + L2s;



% Constantes Résonateur
Lcav =0.16;			% longueur de la cavité
Lcol =0.02;			% longueur du col
Dcav =0.043;			% diametre de la cavité
Dcol =0.02;				% diametre du col
Scav = pi*(Dcav/2)^2;
Scol = pi*(Dcol/2)^2;
RN = Dcol / 2;
RC = Dcav / 2;
RT = d / 2;
L1 = 0.82 * (1 - 1.35*RN/RC + 0.31*(RN/RC)^3) * RN; 
L2 = 0.82 * (1- 0.235 * RN / RT - 1.32*(RN/RT)^2 + 1.54 * (RN/RT)^3 - 0.86*(RN/RT)^4)*RN;
Lcol = Lcol + L1 + L2;

%Base fréquentielle
Fmax=900;
f = 0:0.5:Fmax;
N = length(f);
w = 2*pi*f;
%=================================================================================================================
%Calcul des coefficients de la matrice
reseau = ones(2,2,N);		%matrice de transfert du réseau composé des éléments de 'config.txt'
admittance_resonateur = ones(1,N);


for x=1:1:N
	reseau(:,:,x) = eye(2);
	clc
	fprintf(1,['\b%d  /' num2str(N)],x)	;
	
	for y=1:1:nb_cellule
		if NumSing == y
			matrice_resonateur = resonateur(w(x),Lcavs,Lcols,Dcavs,Dcols,rho,c);		
		else
			matrice_resonateur = resonateur(w(x),Lcav,Lcol,Dcav,Dcol,rho,c);
		end
		
		cellule = guide(w(x),L,d,rho,c) * matrice_resonateur;
		reseau(:,:,x) = reseau(:,:,x)*cellule;
	end
	
	admittance_resonateur(x) = matrice_resonateur(2,1);
end
clc

A = reseau(1,1,:);
B = reseau(1,2,:);
C = reseau(2,1,:);
D = reseau(2,2,:);

%====================================================================
% Affichage dans le terminal
%====================================================================
more off
disp('===============================================================');
if (NumSing>nb_cellule)
	disp("La singularité ne sera pas prise en compte car NumSing>nb_cellule");
end
disp('');
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
disp(['Paramètres Singu']);
disp('---------------------');
disp(['Lcavs = ' num2str(Lcavs) ' m']);
disp(['Lcols = ' num2str(Lcols) ' m']);
disp(['Dcavs = ' num2str(Dcavs) ' m']);
disp(['Dcols = ' num2str(Dcols) ' m']);
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
subplot(3,1,1)
plot(f,cosnGammad,'-b');
hold on
plot(ones(1,N), '-r');
hold on
plot(-ones(1,N), '-r');
axis([0 Fmax -10 10])
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
axis([0 Fmax 0 1.5]);
ylabel('abs(T)');
title('Coefficient de transmission a l entree du reseau');
grid on



%Affichage de l'admittance du résonateur
%-----------------------------------------------
subplot(4,1,4);
semilogy(f,abs(admittance_resonateur(1,:)));
ylabel('abs(Yreso)');
xlabel('frequence en Hz');
title('admittance du resonateur de Helmholtz');
grid on






%Affichage de l'absorption
%-----------------------------------------------
A=1-abs(T).^2 -abs(R).^2;
figure(3)
subplot(2,1,1)
plot(f,abs(T),'b', 'Linewidth',3);
hold on
plot(f,abs(R),'r', 'Linewidth',3);
hold on
plot(f,A,'k');
legend('Transmission','Reflexion','Admittance');
xlim([0 Fmax]);


admittance_singu = ones(1,N);
for x=1:1:N
	matrice_resonateur = resonateur(w(x),Lcavs,Lcols,Dcavs,Dcols,rho,c);		
	admittance_singu(x) = matrice_resonateur(2,1);

end

subplot(2,1,2);
semilogy(f,abs(admittance_singu(1,:)));
grid on
hold on
semilogy(f,abs(admittance_resonateur(1,:)),'-r');
legend('admi singu','admi reso');
xlim([0 Fmax]);
