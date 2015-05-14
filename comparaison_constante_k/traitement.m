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

%[Zc k] = pertes(d,w,rho,c);
%pt 12 du N singu
coef = 28.0498;
coef = 0.074079;

theo = exp(-coef*axe_distance);
axe_distance = axe_distance * 100;   % Passage en cm pour coller au mesures

%======================================================================
%Affichage
figure(1)
plot(defect_high(:,1), defect_high_norm, 'o-r', 'Linewidth',3)
hold on
plot(defect_low(:,1), defect_low_norm,'Linewidth',3);
legend('defect with high lvl','defect with low lvl');
xlabel('distance source recepteur en m');
ylabel('Pression RMS mesuree');
%print -dpng comparaison_lvl_defaut.png

figure(2)
plot(no_defect_high(:,1), no_defect_high_norm, 'o-r','Linewidth',3)
hold on
plot(no_defect_low(:,1), no_defect_low_norm,'Linewidth',3);
legend('no_defect with high lvl','no_defect with low lvl');
xlabel('distance source recepteur en m');
ylabel('Pression RMS mesuree');
%print -dpng comparaison_lvl_sans_defaut.png

figure(3)
plot(defect_low(:,1), defect_low_norm,'o-r','Linewidth',3);
hold on
plot(no_defect_low(:,1), no_defect_low_norm,'Linewidth',3);
legend('defect','no defect');
%print -dpng comparaison_decroissance_lin.png
xlabel('distance source recepteur en m');
ylabel('Pression RMS mesuree');


figure(4)
semilogy(defect_low(:,1), defect_low_norm,'o-r','Linewidth',3);
hold on
semilogy(no_defect_low(:,1), no_defect_low_norm,'Linewidth',3);
legend('defect','no defect');
xlabel('distance source recepteur en m');
ylabel('Pression RMS mesuree');
%print -dpng comparaison_decroissance_log.png

figure(5)
%semilogy(defect_low(:,1), defect_low_norm,'o-r','Linewidth',3);
%hold on
semilogy(no_defect_low(:,1), no_defect_low_norm,'r','Linewidth',3);
hold on
semilogy(axe_distance,theo, 'Linewidth',3);
