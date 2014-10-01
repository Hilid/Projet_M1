function matrice_trans_res = resonateur(k,Lcav,Lcol,Scav,Scol,rho,c) 
% k est un scalaire (calcul pour 1 fréquence)
% Lcav est la longueur de la cavité
% Lcol est la longueur du col
% Scav est la section de la cavité
% Scol est la section du col
%rho est la masse volumique du milieu et c sa célérité

%Fonction qui donne la matrice de transfert d'un résonateur dont l'impédance du fond est infinie

matrice_resonateur = guide(k,Lcol,Scol,rho,c) * guide(k,Lcav,Scav,rho,c);

Zresonateur = matrice_resonateur(1,1)/matrice_resonateur(2,1);	%application de la condition limite : pression max, vitesse nulle.

matrice_trans_res = [ 1, 0 ; 1/Zresonateur, 1];   

endfunction

