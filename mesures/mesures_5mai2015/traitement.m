% Script de traitement des mesures du 05 mai  

%11 résonateurs (5 de part et d'autre du défaut) 
%source HP à -10cm.
%Fréquence d'excitation f=381 Hz
% Le zero se trouve sur la position du défaut
%----------------------------------------------------------------------

clear all; close all;clc;


defect_p=load('with_defect__positive.txt'); 
defect_n=load('with_defect__negative.txt');
defect = [flipud(defect_n); defect_p];            % reconstruction des données

nodefect_p=load('no_defect__positive.txt');
nodefect_n=load('no_defect__negative.txt');
nodefect = [flipud(nodefect_n); nodefect_p];


distance_defect = defect(:,1);
Amp_defect = defect(:,2);

distance_nodefect = nodefect(:,1);
Amp_nodefect = nodefect(:,2);

%-----------------------------------------------------------------------
% Affichage des courbes

figure(1)
plot(distance_defect,Amp_defect./max(Amp_defect),'*-');
hold on
plot(distance_nodefect,Amp_nodefect./max(Amp_nodefect),'r*-');
legend('Avec un defaut','Sans defaut');
title('Champ de pression normalise dans le tube avec et sans defaut');


figure(2)
plot(distance_defect,Amp_defect,'*-');
hold on
plot(distance_nodefect,Amp_nodefect,'r*-');
legend('Avec un defaut','Sans defaut');
title('Champ de pression non normalise dans le tube avec et sans defaut');
