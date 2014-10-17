function matrice_trans_res = resonateur(k,Lcav,Lcol,Dcav,dcol,rho,c) 
% k est un scalaire (calcul pour 1 fréquence)
% Lcav est la longueur de la cavité
% Lcol est la longueur du col
% Dcav est le diametre de la cavité
% dcol est le diametre du col du col
%rho est la masse volumique du milieu et c sa célérité

%Fonction qui donne la matrice de transfert d'un résonateur dont l'impédance du fond est infinie

Scav = pi*(Dcav/2)^2;
Scol = pi*(dcol/2)^2;

%freshelmholtz = c/(2*pi)*sqrt(Scol/(Scav*Lcav*Lcol))    affichage de la fréquence de résonnance

matrice_resonateur = guide(k,Lcol,dcol,rho,c) * guide(k,Lcav,Dcav,rho,c);

Zresonateur = matrice_resonateur(1,1)/matrice_resonateur(2,1);	%application de la condition limite : pression max, vitesse nulle.

matrice_trans_res = [ 1, 0 ; 1/Zresonateur, 1]; 

  

endfunction

