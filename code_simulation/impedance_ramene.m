function a = impedance_ramene(x,L,rho,c,S,Zt)

Zc = rho*c/S;
Zramenee = @(x) Zc *  (Zt + j*Zc*tan(x*L))/(Zc + j * Zt * tan(x*L))  ;   


a = Zramenee(x);
endfunction


%Formule de l'impédance ramenée, utilisée souvent....   Zt => l'impédance a ramener
