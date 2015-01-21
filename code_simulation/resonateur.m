function matrice_trans_res = resonateur(w,Lcav,Lcol,Dcav,dcol,rho,c,d_tuyau) 
% k est un scalaire (calcul pour 1 fréquence)
% Lcav est la longueur de la cavité
% Lcol est la longueur du col
% Dcav est le diametre de la cavité
% dcol est le diametre du col du col
%rho est la masse volumique du milieu et c sa célérité
%d_tuyau le diametre du tuyau sur lequel est fixé le résonateur



%Fonction qui donne la matrice de transfert d'un résonateur dont l'impédance du fond est infinie

Scav = pi*(Dcav/2)^2;
Scol = pi*(dcol/2)^2;

RN = dcol / 2;
RC = Dcav / 2;
RT = d_tuyau / 2;

%Correction de longueur du col (prise dans [A1] appendice B)
L1 = 0.82 * (1 - 1.35*RN/RC + 0.31*(RN/RC)^3) * RN; 
L2 = 0.82 * (1- 0.235 * RN / RT - 1.32*(RN/RT)^2 + 1.54 * (RN/RT)^3 - 0.86*(RN/RT)^4)*RN;
Lcol = Lcol + L1 + L2;


matrice_resonateur = guide(w,Lcol,dcol,rho,c) * guide(w,Lcav,Dcav,rho,c);    %modélisation du résonateur par 2 matrices de guides
Yresonateur = matrice_resonateur(2,1)/matrice_resonateur(1,1);				 %application de la condition limite : pression max, vitesse nulle au fond du résonateur.
matrice_trans_res = [ 1, 0 ; Yresonateur, 1];      						     %matrice du résonateur dans le réseau

  

endfunction

