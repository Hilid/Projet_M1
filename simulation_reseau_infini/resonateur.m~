function matrice_trans_res = resonateur(w,Lcav,Lcol,Dcav,Dcol,rho,c,d_tuyau) 


matrice_resonateur = guide(w,Lcol,Dcol,rho,c) * guide(w,Lcav,Dcav,rho,c);    %modélisation du résonateur par 2 matrices de guides
Yresonateur = matrice_resonateur(2,1)/matrice_resonateur(1,1);				 %application de la condition limite : pression max, vitesse nulle au fond du résonateur.
matrice_trans_res = [ 1, 0 ; Yresonateur, 1];      						     %matrice du résonateur dans le réseau

  

endfunction

