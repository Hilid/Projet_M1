function matrice_trans_guide = guide(w,L,d,rho,c)
R = d/2; %rayon du guide
S = pi*(R)^2;     %section du guide

%Cas sans perte
%--------------
Zc = rho*c/S;  %impedance caracteristique du guide
k = w/ c;

%Cas avec pertes (page 4 [A1])
%-----------------------------
%[Zc k] = pertes(d,w,rho,c);



matrice_trans_guide =  [cos(k*L)  Zc*j*sin(k*L) ; (1/Zc)*j*sin(k*L) cos(k*L)];   

endfunction
