clear all
close all
clc
graphics_toolkit('gnuplot')          %affichage gnuplot
more off
%hold on

%===============================================================================================================
%Constantes physiques
%--------------------
%Fréquence imposée
f= 371; 								%866 : Bragg  1000 : R=1
w=2*pi*f;

%Périodisation
nb_cellule =5;

%Nombre de points de visualisation pour chaque guide
Ndiv = 50;   							%les portions de guide sont divisées en 5 sous-guides

%Singularité sur la longueur de cavité du résonateur numéro NumSing
NumSing= 3;

%Constantes milieu
c = 343;
rho = 1.177;    		   %a  300°K

% Constantes guide
L = 0.1; 					% longueur du guide
d = 0.05; 					% diametre du guide
L_div = L/Ndiv; 						%longueur d'un sous-guide

%Constantes Résonateur
Lcav =0.16;					% longueur de la cavité
Lcol =0.02;					% longueur du col
Dcav =0.043;				% diametre de la cavité
Dcol =0.02;					% diametre du col

Scav = pi*(Dcav/2)^2;
Scol = pi*(Dcol/2)^2;

RN = Dcol / 2;
RC = Dcav / 2;
RT = d / 2;

%Correction de longueur du col (prise dans [A1] appendice B)
L1 = 0.82 * (1 - 1.35*RN/RC + 0.31*(RN/RC)^3) * RN; 
L2 = 0.82 * (1- 0.235 * RN / RT - 1.32*(RN/RT)^2 + 1.54 * (RN/RT)^3 - 0.86*(RN/RT)^4)*RN;
Lcol = Lcol + L1 + L2;


Lcavs =0.08;  					
Lcols =0.02;					
Dcavs =0.043;				
Dcols =0.02;					

Scavs = pi*(Dcavs/2)^2;
Scols = pi*(Dcols/2)^2;

RNs = Dcols / 2;
RCs = Dcavs / 2;
RTs = d / 2;

%Correction de longueur du col (prise dans [A1] appendice B)
L1s = 0.82 * (1 - 1.35*RNs/RCs + 0.31*(RNs/RCs)^3) * RNs; 
L2s = 0.82 * (1- 0.235 * RNs / RTs - 1.32*(RNs/RTs)^2 + 1.54 * (RNs/RTs)^3 - 0.86*(RNs/RTs)^4)*RNs;
Lcols = Lcols + L1s + L2s;


%Matrice pression-vitesse imposée à la sortie
PV_out = [1 ; (d/2)^2*pi/(rho*c)]; 						 %sortie anéchoïque 


%Nombre de matrice de transfert et donc de point de visualisation total
NbPtTot = nb_cellule*(Ndiv +1); % +1 pour chaque résonateur

%=================================================================================================================
%Calcul des coefficients de la matrice
%--------------------------------------
loc=NbPtTot;					%localisation du point de visualisation dans le guide. On part de la sortie.
PV=zeros(2,NbPtTot+1);			%initialisation de la matrice contenant pression et vitesse (en colonne) pour chaque point (NbPtTot lignes + 1 pour la sortie)
PV(:,loc+1) = PV_out; 			%matrice pression-vitesse en UN SEUL point du réseau, pour les N fréquences

for y=1:1:nb_cellule

	if (NumSing==(nb_cellule+1)-y)
		PV(:,loc)=resonateur(w,Lcavs,Lcols,Dcavs,Dcols,rho,c)*PV(:,loc+1);
	else
		PV(:,loc)=resonateur(w,Lcav,Lcol,Dcav,Dcol,rho,c)*PV(:,loc+1);
	end

	loc=loc-1;					%on avance d'un point (de la sortie vers l'entrée)
	for ndiv=1:Ndiv
		PV(:,loc) = guide(w,L_div,d,rho,c)*PV(:,loc+1);
		loc=loc-1;  
	end
	
end

%====================================================================
% Affichage dans le terminal
%----------------------------
disp('===============================================================');
if (NumSing>nb_cellule|| NumSing<0)
	disp("La singularité ne sera pas prise en compte car NumSing>nb_cellule");
end
disp(['Paramètres globaux']);
disp('------------------');
disp(['frequence imposée = ' num2str(f) ' hz']);
disp(['célérité: ',num2str(c) ' m.s^-1']);
disp(['masse volumique: ',num2str(rho) ' kg.m^-3']);
disp('');
disp(['Paramètres résonateur']);
disp('---------------------');
disp(['Lcav = ' num2str(Lcav) ' m']);
disp(['Lcavs = ' num2str(Lcavs) ' m']);
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
disp(['Paramètres réseau']);
disp(['nb_cellule = ' num2str(nb_cellule) ]);
disp(['nombre de division par guide = ' num2str(Ndiv)]);
disp('===============================================================');
disp('');

%================================================================================
%Affichage de la pression
%-------------------------
figure(1)
plot(abs(real((PV(1,:)))),'o-')
ylabel('Pression')
%ylim([0 10^10]);
title(['Pression dans le guide pour f=',num2str(f)]);
hold on
Pmin=min(abs(real(PV(1,:))));
Pmax=max(abs(real(PV(1,:))));
hold on
for n=1:nb_cellule
	line((Ndiv+1)*n+1,[Pmin Pmax]); %ligne avant résonateur
	line((Ndiv+1)*n,[Pmin Pmax]);	%ligne après résonateur
end
hold off

%Affichage de la vitesse
%-----------------------
%~ figure(5)
%~ plot(real(PV(2,:)),'o-')
%~ ylabel('Vitesse')
%~ title(['Vitesse dans le guide pour f=',num2str(f)]);
%~ Pmin=min(real(PV(2,:)));
%~ Pmax=max(real(PV(2,:)));
%~ hold on
%~ for n=1:nb_cellule
	%~ line((Ndiv+1)*n+1,[Pmin Pmax]); %ligne avant résonateur
	%~ line((Ndiv+1)*n,[Pmin Pmax]);	%ligne après résonateur
%~ end
%~ hold off
