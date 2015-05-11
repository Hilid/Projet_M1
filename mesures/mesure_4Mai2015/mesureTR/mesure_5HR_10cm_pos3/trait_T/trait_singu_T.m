clear all
close all
clc
graphics_toolkit('gnuplot')          %affichage gnuplot

%===============================================================================================================
nb_cellule =5;

%Constantes physiques
c = 340;
rho = 1.177;       %a  300°K

% Constantes guide
L = 0.10; 				% longueur du guide
d = 0.05; 				% diametre du guide

% Constantes Singularité
NumSing = 3;

Lcavs =0.095;			% seul truc qui change
Lcols =0.02;			
Dcavs =0.043;			
Dcols =0.02;				
Scavs = pi*(Dcavs/2)^2;
Scols = pi*(Dcols/2)^2;
RNs = Dcols / 2;
RCs = Dcavs / 2;
RTs = d / 2;
L1s = 0.82 * (1 - 1.35*RNs/RCs + 0.31*(RNs/RCs)^3) * RNs; 
L2s = 0.82 * (1- 0.235 * RNs / RTs - 1.32*(RNs/RTs)^2 + 1.54 * (RNs/RTs)^3 - 0.86*(RNs/RTs)^4)*RNs;
Lcols = Lcols + L1s + L2s;



% Constantes Résonateur
Lcav =0.16;			% longueur de la cavité
Lcol =0.02;			% longueur du col
Dcav =0.043;			% diametre de la cavité
Dcol =0.02;				% diametre du col
Scav = pi*(Dcav/2)^2;
Scol = pi*(Dcol/2)^2;
RN = Dcol / 2;
RC = Dcav / 2;
RT = d / 2;
L1 = 0.82 * (1 - 1.35*RN/RC + 0.31*(RN/RC)^3) * RN; 
L2 = 0.82 * (1- 0.235 * RN / RT - 1.32*(RN/RT)^2 + 1.54 * (RN/RT)^3 - 0.86*(RN/RT)^4)*RN;
Lcol = Lcol + L1 + L2;

%Base fréquentielle
Fmax=900;
f = 150:1:Fmax;
N = length(f);
w = 2*pi*f;
%=================================================================================================================
%Calcul des coefficients de la matrice
reseau = ones(2,2,N);		%matrice de transfert du réseau composé des éléments de 'config.txt'
admittance_resonateur = ones(1,N);


for x=1:1:N
	reseau(:,:,x) = eye(2);	
	for y=1:1:nb_cellule
		if NumSing == y
			matrice_resonateur = resonateur(w(x),Lcavs,Lcols,Dcavs,Dcols,rho,c);		
		else
			matrice_resonateur = resonateur(w(x),Lcav,Lcol,Dcav,Dcol,rho,c);
		end
		
		cellule = guide(w(x),L,d,rho,c) * matrice_resonateur;
		reseau(:,:,x) = reseau(:,:,x)*cellule;
	end
	
	admittance_resonateur(x) = matrice_resonateur(2,1);
end
clc

A = reseau(1,1,:);
B = reseau(1,2,:);
C = reseau(2,1,:);
D = reseau(2,2,:);

Zr = sqrt(B./C);


S = d^2/4*pi;
Zc = rho*c/S;

Zc_perte = ones(1,1,N);

for x=1:1:N
		[Zc_perte(1,1,x) b] = pertes(d,w(x),rho,c);
end
Zc = Zc_perte;


R = (A + B./Zc - C.*Zc - D)./(A + C.*Zc + B./Zc + D);


T = 2./(A + C.*Zc + B./Zc + D);


%%%%%%%%
%plot
%%%%%%%%
figure(1)
plot(f,abs(T(1,1,:)));
axis([0 Fmax 0 1.5]);
ylabel('|T|');
xlabel('Frequences (Hz)')
%title('Coefficient de transmission a l entree du reseau');
grid on
axis([min(f) max(f) 0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DONNEES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% geometrical parameters of the sensor and of the tube
d1  = 18*0.001;        %diamétre cavité fermee
s1  = pi*(d1/2)^2;     %section cavité fermee
d2  = 16*0.001;        %diamétre cavite ouverte
s2  = pi*(d2/2)^2;     %section cavité ouverte
l1  = 21.4*0.001;      %longueur de la cavite fermee
l12 = 9.5*0.001;       %demi longueur de la cavite fermee 
l2  = 13*0.001;        %longueur de la cavite ouverte 
l22 = 6.5*0.001;       %demi longueur de la cavite ouverte

R   = (3.3/2)*10^(-2);     %radius of the tube 
S   = pi*R^2;          %cross section area

%% Air properties
T_0    = 20;                                     %[°C]
P_0    = 101320;                                 % ambiant pressure [Pa]
mmol   = 28.96e-3;                               %[kg.mol^-1]
gamma  = 1.4;                     % ratio of specific heats, le rapport des chaleurs specifiques de l'air
rho    = P_0*mmol/(8.3123*(273.15+T_0));         %[kg.m^-3]
c      = sqrt(gamma*8.3143*(273.15+T_0)/mmol);   %[m.s^-1]
Cp     = 1006;                     % heat capacity under constant pressure for air
Pr     = 0.71;                       % at atmosperic pressure, C_p*?/K
mu     = 18.4e-6;                  % Dynamic viscosity
kappa  = 2.62*10^(-2);              % conductivité thermique [W/m/K]
nu     = 15.6*10^(-6);

%% Characteristic Impedances and quantities
load ../H13_cal.txt
H13_cal = H13_cal;
ReH13_cal  = H13_cal(:,2);
ImH13_cal  = H13_cal(:,3);
H13_cal    = ReH13_cal + j*(ImH13_cal); 


load ../H13.txt
ReH13    = H13(:,2);
ImH13    = H13(:,3);
H13      = ReH13 + j*ImH13;   %Transfer function between the first and the third microphone

 
load ../Impedance.txt
f     = Impedance(:,1);
ReImp = Impedance(:,2);
ImImp = Impedance(:,3);
Z     = ReImp+j*ImImp;    %Entrance impedance

omega  = 2*pi*f;

k      = omega./c;
Zt     = (rho*c/S);
Zc1    = rho*c./s1;       %characteristic impedance de la cavity fermee
Zc2    = rho*c./s2;       %characteristic impedance de la cavity ouverte

%K      = -j.*1./Zc1.*(1./H12_cal).*sin(k.*l1).*cos(k.*l22)./cos(k.*l12)./cos(k.*l2);
%beta   = j.*Zc2.*tan(k.*l22);
delta  = j.*tan(k.*l2)./Zc2;

%Zn     = (H12./K-beta)./(1-delta.*H12./K);  % Entrance impedance
K_T    =  H13_cal;
l=0.05;
H13=H13.*exp(j*k*l);
%-H13_cal.*Zc2./Zc1.*sin(k.*l1)./(sin(k.*l2).*cos(k.*l12));
Z_T    =  2*H13./delta./K_T.*(1+delta.*Z.*Zt)./Zt;  %transfer impedance


figure(200)
subplot(211)
plot(omega./2./pi,20*log10(abs(Z)),'k')
set(gca,'Fontsize',15)
xlabel('freq(Hz)','Fontsize',18)
ylabel('|Z_{in}| (dB)','Fontsize',18)
hold on
subplot(212)
plot(omega./2./pi,(angle(Z)),'k')
hold on
set(gca,'Fontsize',15)
xlabel('freq(Hz)','Fontsize',18)
ylabel('angle(Z_{in})','Fontsize',18)


save Z_T 'Z_T'

T = 2*Z_T.*Zt./(Z.*Zt+Zt);
Re = (Z.*Zt-Zt)./(Z.*Zt+Zt);
Ab = 1-abs(T).^2-abs(Re).^2;

% phi  = -unwrap(phase(T));
% f1   = 300:.8:800;
% phi1 = spline(omega./2./pi,phi,f1);
% phi1 = smooth(phi1,10);
% omega1 = 2.*pi.*f1;
% inv_group = diff(phi1)./diff(omega1');


%%%%%%%%%%%%
%plot
%%%%%%%%%%
figure(1)
hold on
plot(omega./2./pi,(abs(T)),'k')%,omega./2./pi,20*log10(Ab),'r')





