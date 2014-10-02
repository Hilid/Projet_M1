clear all
close all
clc
%graphics_toolkit('gnuplot')          %affichage gnuplot

%===================================================================
%Chargement du fichier de configuration
%--------------------------------------
donnees = load('config.txt');
nb_element = size(donnees)(1);


%Constantes 
%----------
c = 340;
rho = 1.2;

%Base fréquentielle
%------------------
Fmax=2000;
f = 0:0.1:Fmax;
w = f .* (2*pi);
k = w / c;
N = length(f);

%====================================================================

%Calcul des coefficients de la matrice
%--------------------------------------
res = ones(1,N);			%matrice de transfert du résonateur
reseau = ones(2,2,N);		%matrice de transfert du réseau composé des éléments de 'config.txt'


for x=1:1:N
	reseau(:,:,x) = eye(2);	%initialisation de la matrice
end


for x=1:1:N
	y = k(x);
	for d=1:1:nb_element
		if ((donnees(d,1)==1))
			reseau(:,:,x) = reseau(:,:,x) * guide(y,donnees(d,2), donnees(d,3),rho,c);
		elseif ((donnees(d,1)==2))
		     reseau(:,:,x) =reseau(:,:,x) * resonateur(y,donnees(d,2),donnees(d,3),donnees(d,4),donnees(d,5),rho,c); 
		     admittance_resonateur(x) = resonateur(y,donnees(2,2),donnees(2,3),donnees(2,4),donnees(2,5),rho,c)(2,1);
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


%Affichage de l'équation de dispersion 2cos(Gamma d) = T11 + T12
%--------------------------------------------------------------
cosGammad = ((reseau(1,1,:) + reseau(2,2,:)) /2);

figure(1)
plot(f,cosGammad,'-b');
hold on
plot(ones(1,N), '-r');
hold on
plot(-ones(1,N), '-r');
axis([0 1000 -10 10])
legend('cos( \Gamma d )','y=1','y=-1');
xlabel('frequence en Hz');
ylabel('cos( \Gamma d)');
title("Representation graphique de l'equation de dispersion")



%Affichage de l'admittance du résonateur
%--------------------------------------------------------------
figure(2)
subplot(3,1,1)
plot(f,log(abs(admittance_resonateur)));
%xlabel('frequence en Hz');
ylabel('admittance');
title('admittance du resonateur de Helmholtz');

%Affichage de l'impédance caractéristique du réseau Zc
%-----------------------------------------------------
Zc = sqrt(reseau(1,2,:)./reseau(2,1,:));

hold on
plot(f,log(abs(Zc)), 'r');
%xlabel('frequence en Hz');
ylabel('abs(Zcreseau)');
title('Impedance caracteristique du reseau');
grid minor on

%Affichage du coefficient de reflexion
%--------------------------------------------------------------
R = (reseau(1,1,:) + reseau(1,2,:)./Zc - reseau(2,1,:).*Zc- reseau(2,2,:))./(reseau(1,1,:) + reseau(2,1,:).*Zc + reseau(1,2,:)./Zc + reseau(2,2,:));

figure(2)
subplot(3,1,2);
plot(f,abs(R));
%xlabel('frequence en Hz');
ylabel('abs(R)');
title('Coefficient de reflexion a l entree du reseau')
grid minor on

%Affichage du coefficient de transmission
%--------------------------------------------------------------
T = 2./(reseau(1,1,:)  + reseau(2,1,:).*Zc + reseau(1,2,:)./Zc + reseau(2,2,:));

figure(2)
subplot(3,1,3);
plot(f,abs(T));
%xlabel('frequence en Hz');
ylabel('abs(T)');
title('Coefficient de transmission a l entree du reseau');
grid minor on
