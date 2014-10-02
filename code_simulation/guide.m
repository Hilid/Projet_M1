function matrice_trans_guide = guide(k,L,d,rho,c)
% k est un scalaire (calcul pour 1 fréquence)
% L est la longueur du guide
% d est le diametre du guide
%rho est la masse volumique du milieu et c sa célérité

%Fonction qui donne la matrice de transfert d'un guide d'onde (ou d'un fil)
S = pi*(d/2)^2; %section du guide
Zc = rho*c/S;  %impedance caracteristique du guide

matrice_trans_guide =  [cos(k*L)  Zc*j*sin(k*L) ; (1/Zc)*j*sin(k*L) cos(k*L)];   

endfunction
