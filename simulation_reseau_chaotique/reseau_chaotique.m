clear all
close all
clc
graphics_toolkit('gnuplot')          %affichage gnuplot

%===============================================================================================================
%Chargement du fichier de configuration
%--------------------------------------
donnees = load('config.txt');
nb_element = size(donnees)(1);


%Constantes 
%----------
c = 340;
rho = 1.177;       %a  300°K

%Base fréquentielle
%------------------
Fmax=2000;
f = 0:1:Fmax;
N = length(f);

%=================================================================================================================

%Calcul des coefficients de la matrice
%--------------------------------------
res = ones(1,N);			%matrice de transfert du résonateur
reseau = ones(2,2,N);		%matrice de transfert du réseau composé des éléments de 'config.txt'


for x=1:1:N
	reseau(:,:,x) = eye(2);	%initialisation de la matrice par une matrice diago de 1
end


for x=1:1:N
	w = 2*pi* x / N * Fmax ;
	for d=1:1:nb_element
		if ((donnees(d,1)==1))
			reseau(:,:,x) = reseau(:,:,x) * guide(w,donnees(d,2), donnees(d,3),rho,c);
			d_tuyau = donnees(d,3);                                  %sauvegarde du dernier diametre donné pour avoir la section actuelle du guide
		elseif ((donnees(d,1)==2))
			 cellule_reso = resonateur(w,donnees(d,2),donnees(d,3),donnees(d,4),donnees(d,5),rho,c,d_tuyau);
		     reseau(:,:,x) =reseau(:,:,x) * cellule_reso; 
		     admittance_resonateur(x) = cellule_reso(2,1);
		else
			 disp('erreur: element non reconnu');
			 return
		end
	end
end

%====================================================================
% Affichage des paramètres globaux
%----------------------------------

disp('===============================================================');
disp(['Paramètres globaux']);
disp('------------------');
disp(['célérité: ',num2str(c)]);
disp(['masse volumique: ',num2str(rho)]);
disp(['nombre d éléments du réseau : ',num2str(nb_element)]);
disp('===============================================================');

%Affichage de l'impédance caractéristique du réseau Zc
%-----------------------------------------------------
Zc = sqrt(reseau(1,2,:)./reseau(2,1,:));

figure(2)
subplot(4,1,1);
plot(f,log(abs(Zc)), 'r');
%xlabel('frequence en Hz');
ylabel('abs(Zcreseau)');
title('Impedance caracteristique du reseau');
grid minor on

%Affichage du coefficient de reflexion
%--------------------------------------------------------------
R = (reseau(1,1,:) + reseau(1,2,:)./Zc - reseau(2,1,:).*Zc- reseau(2,2,:))./(reseau(1,1,:) + reseau(2,1,:).*Zc + reseau(1,2,:)./Zc + reseau(2,2,:));

figure(2)
subplot(4,1,2);
plot(f,abs(R));
%xlabel('frequence en Hz');
ylabel('abs(R)');
title('Coefficient de reflexion a l entree du reseau')
grid minor on

%Affichage du coefficient de transmission
%--------------------------------------------------------------
T = 2./(reseau(1,1,:)  + reseau(2,1,:).*Zc + reseau(1,2,:)./Zc + reseau(2,2,:));

figure(2)
subplot(4,1,3);
plot(f,abs(T));
%xlabel('frequence en Hz');
ylabel('abs(T)');
title('Coefficient de transmission a l entree du reseau');
grid minor on



%Affichage de l'admittance du dernier résonateur
%-----------------------------------------------
figure(2)
subplot(4,1,4);
plot(f,log(abs(admittance_resonateur)));
ylabel('abs(Yreso)');
xlabel('frequence en Hz');
title('admittance du resonateur de Helmholtz');
grid minor on

