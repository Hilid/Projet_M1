function matrice_trans_guide = guide(w,L,d,rho,c)
% k est un scalaire (calcul pour 1 fréquence)
% L est la longueur du guide
% d est le diametre du guide
%rho est la masse volumique du milieu et c sa célérité

%Fonction qui donne la matrice de transfert d'un guide d'onde (ou d'un fil)

S = pi*(d/2)^2;     %section du guide
R = d/2;            %diametre du guide

%Cas sans perte
%--------------
Zc = rho*c/S;  %impedance caracteristique du guide
k = w/ c;

%Cas avec pertes (page 4 [A1])
%-----------------------------
%~ Pr = 0.708;                                        % nombre de Prandtl a la pression atmosphérique           a 300°K
%~ mu = 1.85 * 10^-5;                                       % viscosité de l'air                                a 300°K
%~ gamma = 1.4;                                   % heat capacity ratio of air
%~ 
%~ ksi = sqrt(Pr);
%~ delta = sqrt(2 * mu / (rho * w));             % viscous boundary layer thickness
%~ s = R /delta;
%~ beta = (1-i)/sqrt(2);
%~ 
%~ k = w / c * (1 + beta /s * (1 + (gamma - 1)/ksi));
%~ Zc = rho * c / S * (1 + beta /s * (1 - (gamma - 1)/ksi));


matrice_trans_guide =  [cos(k*L)  Zc*j*sin(k*L) ; (1/Zc)*j*sin(k*L) cos(k*L)];   

endfunction
