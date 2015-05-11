function [Zc k] = pertes(d,w,rho,c);
R = d/2; 									% Rayon du tube considéré
S = pi*(R)^2;    


Pr = 0.708;                                        % nombre de Prandtl a la pression atmosphérique           a 300°K
mu = 1.85 * 10^-5;                                 % viscosité de l'air		                                a 300°K
gamma = 1.4;                               	      % heat capacity ratio of air

ksi = sqrt(Pr);
delta = sqrt(2 * mu / (rho * w));            	  % viscous boundary layer thickness
s = R /delta;
beta = (1-i)/sqrt(2);

Zc = rho * c / S * (1 + beta /s * (1 - (gamma - 1)/ksi));
k = w / c * (1 + beta /s * (1 + (gamma - 1)/ksi));

endfunction
