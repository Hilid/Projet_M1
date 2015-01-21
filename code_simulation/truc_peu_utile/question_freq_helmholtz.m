% Pourquoi l'approximation basses fréquences donne-t-elle une fréquence de résonance aussi éloignée de la solution analytique ?
clear all; close all;

D = 0.05;	%diamètre cavité
L = 0.165;	%hauteur cavité
S = pi*(D/2)^2;	%section cavité
V = L*S;	%volume cavité

d = 0.02; %diamètre goulot
l = 0.02; %longueur du col
s = pi*(d/2)^2; %section du col

c=340;
rho = 1.2;
Zc1 = rho*c/S;	%impédance cavite
Zc2 = rho*c/s;	%impedance volume

fmax = 1600;
f=0:1:fmax;
k= 2*pi*f/c;

% Formule avec approx kl<<1

fres = c*sqrt(s/(V*l))/(2*pi);
disp(["avec la formule, fres = " num2str(fres)]);

% Résolution numérique de la formule analytique

a = Zc1-Zc2*tan(k*l).*tan(k*L);
b= Zc2*j*tan(k*L)+j*Zc1*tan(k*l);

Z = Zc2*(a./b); 

plot(f,zeros(size(f)));
hold on
plot(f,a);
axis([0 fmax -2*10^6 2*10^6]);
title("Dénominateur de l'impédance -> resonance quand =0");

%figure
%semilogy(f,abs(Z));

figure 
semilogy(f(2:end),abs(1./Z(2:end)));
title('Admittance du résonateur');
[valeur fres]=max(abs(1./Z));
disp(["en calculant l'impedance ramenée, fres = " num2str(fres)]);



