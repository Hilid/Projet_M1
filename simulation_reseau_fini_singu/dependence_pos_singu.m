clear all
close all
clc
graphics_toolkit('gnuplot')          %affichage gnuplot

%===============================================================================================================
nb_cellule =5;

%Constantes physiques
c = 340;
rho = 1.177;       %a  300°K

% Constantes guide
L = 0.10; 				% longueur du guide
d = 0.05; 				% diametre du guide

% Constantes Singularité
Lcavs =0.065;			% seul truc qui change
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
Lcav =0.145;			% longueur de la cavité
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
Fmax=600;
f = 100:0.5:Fmax;
N = length(f);
w = 2*pi*f;
%=================================================================================================================
%Calcul des coefficients de la matrice
reseau = ones(2,2,N);		%matrice de transfert du réseau composé des éléments de 'config.txt'
admittance_resonateur = ones(1,N);

for pos_singu=1:nb_cellule
	for x=1:1:N
		reseau(:,:,x) = eye(2);
		matrice_resonateur_singulier = resonateur(w(x),Lcavs,Lcols,Dcavs,Dcols,rho,c);		
		matrice_resonateur = resonateur(w(x),Lcav,Lcol,Dcav,Dcol,rho,c);
		matrice_guide = guide(w(x),L,d,rho,c);
		reseau(:,:,x) = (matrice_guide * matrice_resonateur)^(pos_singu-1) * (matrice_guide * matrice_resonateur_singulier) * ((matrice_guide * matrice_resonateur)^(nb_cellule-pos_singu));
		admittance_resonateur(x) = matrice_resonateur_singulier(2,1);
	end

	A = reseau(1,1,:);
	B = reseau(1,2,:);
	C = reseau(2,1,:);
	D = reseau(2,2,:);
	
	%Affichage du coefficient de transmission
	%--------------------------------------------------------------
	S = d^2/4*pi;
	Zc = rho*c/S;

	Zc_perte = ones(1,1,N);

	for x=1:1:N
		[Zc_perte(1,1,x) b] = pertes(d,w(x),rho,c);
	end
	Zc = Zc_perte;
	
	T = 2./(A + C.*Zc + B./Zc + D);
	vec_T(pos_singu,:)=T(1,1,:);
end

figure
subplot(211)
plot(f,abs(vec_T(3,:)))
subplot(212)
plot(f,admittance_resonateur)

figure
mesh(f,1:nb_cellule, abs(vec_T));
view(0,90)
