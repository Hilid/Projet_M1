clear all
close all
clc
format long

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
load H13_cal.txt
H13_cal = H13_cal;
ReH13_cal  = H13_cal(:,2);
ImH13_cal  = H13_cal(:,3);
H13_cal    = ReH13_cal + j*(ImH13_cal); 
%H13_cal    = abs(H13_cal).*exp(j.*unwrap(angle(H13_cal)));

% load H12_cal.txt
% ReH12_cal  = H12_cal(:,2);
% ImH12_cal  = H12_cal(:,3);
% H12_cal    = ReH12_cal + j*ImH12_cal; 

load H13.txt
ReH13    = H13(:,2);
ImH13    = H13(:,3);
H13      = ReH13 + j*ImH13;   %Transfer function between the first and the third microphone

% load H12.txt
% ReH12    = H12(:,2);
% ImH12    = H12(:,3);
% H12      = ReH12 + j*ImH12;   %Transfer function between the first and the second microphone

% load P3surQ.txt
% Re      = P3surQ(:,2);
% Im      = P3surQ(:,3);
% P3surQ  = (rho*c/S).*(Re + j*Im);   
% 
load Impedance.txt
f     = Impedance(:,1);
ReImp = Impedance(:,2);
ImImp = Impedance(:,3);
Z     = ReImp+j*ImImp;    %Entrance impedance

omega  = 2*pi*f;
% Tube
% rv_t   = R*sqrt(omega.*rho./mu);                           % le rayon du tube incluant les pertes visqueuses
% rt_t   = R*sqrt(omega.*rho.*Cp./kappa);                    % le rayon du tube incluant les pertes thermiques
% Zv     = j.*omega.*rho./S.*(1+sqrt(2)./rv_t.*(1-j));       % l'impedance de la tube incluant les pertes
% Yt     = j.*omega.*S./rho./c^2.*(1+(gamma-1)*sqrt(2)./rt_t.*(1-j)); %l'admittance
% Gamma  = sqrt(Zv.*Yt);
% kt     = -j*Gamma;
%Zt     = sqrt(Zv./Yt);   % characteristic impedance of the tube
k      = omega./c;
Zt     = (rho*c/S);
Zc1    = rho*c./s1;       %characteristic impedance de la cavity fermee
Zc2    = rho*c./s2;       %characteristic impedance de la cavity ouverte

%K      = -j.*1./Zc1.*(1./H12_cal).*sin(k.*l1).*cos(k.*l22)./cos(k.*l12)./cos(k.*l2);
%beta   = j.*Zc2.*tan(k.*l22);
delta  = j.*tan(k.*l2)./Zc2;

%Zn     = (H12./K-beta)./(1-delta.*H12./K);  % Entrance impedance
K_T    =  H13_cal;
%-H13_cal.*Zc2./Zc1.*sin(k.*l1)./(sin(k.*l2).*cos(k.*l12));
Z_T    =  H13./delta./K_T.*(1+delta.*Z.*Zt)./Zt;  %transfer impedance

figure(100)
subplot(211)
plot(omega./2./pi,20*log10(abs(Z_T)),'k')
set(gca,'Fontsize',15)
xlabel('freq(Hz)','Fontsize',18)
ylabel('|Z_{T}| (dB)','Fontsize',18)
hold on
subplot(212)
plot(omega./2./pi,(angle(Z_T)),'k')
set(gca,'Fontsize',15)
xlabel('freq(Hz)','Fontsize',18)
ylabel('angle(Z_{T})','Fontsize',18)

hold on

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

figure(300)
subplot(211)
plot(omega./2./pi,20*log10(abs(T).^2),'k')%,omega./2./pi,20*log10(Ab),'r')
hold on
xlim([100 1000])
subplot(212)
plot(omega./2./pi,unwrap(angle(T)),'k')
hold on

figure(400)
plot(omega./2./pi,(abs(T).^2),'k') %,omega./2./pi,Ab,'r',omega./2./pi,(abs(Re).^2))
set(gca,'Fontsize',15)
hold on
xlabel('freq(Hz)','Fontsize',18)
ylabel('Transmission','Fontsize',18)
xlim([100 1000])
ylim([0 1.3])
hold on

figure(500)
plot(f,abs(Re))
title('|R|')

%figure(500)
%plot(omega1(1:end-1)./2./pi,inv_group)

% save lc16cmlc26cm 'T' 'Ab' 'Re' 
