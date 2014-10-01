function a = masse(x,M)

sing_masse_corde = @(x) [1 j*x*M ; 0 1];

a = sing_masse_corde(x);

endfunction


% Matrice pour une masse ( de masse M) sur un fil
