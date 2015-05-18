% Script de traitement et d'affichage des mesures de pression du 12 mai
clear all
close all
clc
graphics_toolkit('fltk');


%fréquence d'excitation
F = 287.7;

% Chargement des mesures
defect_low = load('defect_low.txt');
defect_high = load('defect_high.txt');

no_defect_low = load('no_defect_low.txt');
no_defect_high = load('no_defect_high.txt');


% Normalisation
defect_low_norm = defect_low(:,2) ./ max(defect_low(:,2));
defect_high_norm = defect_high(:,2) ./ max(defect_high(:,2));
no_defect_low_norm = no_defect_low(:,2) ./ max(no_defect_low(:,2));
no_defect_high_norm = no_defect_high(:,2) ./ max(no_defect_high(:,2));

%==========================================================================
%Calcul du coefficient de perte correspondant
N = 50;
d_max = 0.8;      %distance max de 80 cm
axe_distance = (1:1:N)./N .*d_max;

coef_defaut =  5.92387; % Coef pour 2 singu

%coef_defaut = 6.77178; % Coef pour une seul singu

theo_avecdefaut = exp(-coef_defaut*axe_distance);

coef_sansdefaut=7.4115;
theo_sansdefaut=exp(-coef_sansdefaut*axe_distance);

axe_distance = axe_distance * 100;   % Passage en cm pour coller au mesures

%======================================================================
%Affichage
figure(1)
plot(defect_high(:,1), defect_high_norm, 'o-r', 'Linewidth',3)
hold on
plot(defect_low(:,1), defect_low_norm,'o-b','Linewidth',3);
legend('defect with high lvl','defect with low lvl');
xlabel('distance source recepteur en m');
ylabel('Amplitude (V)');
%print -dpng comparaison_lvl_defaut.png

figure(2)
plot(no_defect_high(:,1), no_defect_high_norm, 'o-r','Linewidth',3)
hold on
plot(no_defect_low(:,1), no_defect_low_norm,'o-b','Linewidth',3);
legend('no_defect with high lvl','no_defect with low lvl');
xlabel('distance source recepteur en m');
ylabel('Amplitude (V)');
%print -dpng comparaison_lvl_sans_defaut.png

figure(3)
plot(defect_low(:,1), defect_low_norm,'o-r','Linewidth',3);
hold on
plot(no_defect_low(:,1), no_defect_low_norm,'o-b','Linewidth',3);
hold on
legend('Avec defaut','Sans defaut');
xlabel('Distance source recepteur en cm');
ylabel('Amplitude (V)');
print -dpng comparaison_decroissance_lin_theo.png


figure(4)
semilogy(defect_low(:,1), defect_low_norm,'o-r','Linewidth',3);
hold on
semilogy(no_defect_low(:,1), no_defect_low_norm,'o-b','Linewidth',3);
xlabel('Distance source recepteur en cm');
ylabel('Amplitude (V)');

hold on
%semilogy(axe_distance,theo_avecdefaut,'-r', 'Linewidth',2);
hold on
semilogy(axe_distance,theo_sansdefaut,'--b', 'Linewidth',2);
hold on
pos_10cm=ceil(10*N/d_max/100); %positition de 10cm dans axe_distance
semilogy(axe_distance(pos_10cm:end),exp(-coef_sansdefaut*axe_distance(1:end-pos_10cm+1)./100),'r--', 'Linewidth',2)
legend('Mesures avec défaut','Mesures sans défaut','Théorie sans défaut','Théorie sans défaut, partant de x=10cm');



print -dsvg comparaison_decroissance_log_theo.svg
