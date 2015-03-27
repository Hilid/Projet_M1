function f2l(f)

Lcol=0.02;
Lcav =0.16;			% longueur de la cavité
Dcav =0.043;			% diametre de la cavité
Dcol =0.02;				% diametre du col
Scav = pi*(Dcav/2)^2;
Scol = pi*(Dcol/2)^2;
RN = Dcol / 2;
RC = Dcav / 2;
L = 0.10; 				% longueur du guide
d = 0.05; 				% diametre du guide
RT = d / 2;
L1 = 0.82 * (1 - 1.35*RN/RC + 0.31*(RN/RC)^3) * RN; 
L2 = 0.82 * (1- 0.235 * RN / RT - 1.32*(RN/RT)^2 + 1.54 * (RN/RT)^3 - 0.86*(RN/RT)^4)*RN;
Lcol=Lcol+L1+L2;
c=343;

Lcav=  (c^2 * Scol)/ (4*pi*pi*f^2*Scav*Lcol); %longueur equivalente = avec correction

w=2*pi*f;
L=(0.5:0.001:1.5)*Lcav;
N=length(L);


for x=1:N
	matrice_resonateur = resonateur(w,L(x),Lcol,Dcav,Dcol,1.2,c);		
	admittance_singu(x) = matrice_resonateur(2,1);
end


%~ figure
%~ plot(L,(abs(admittance_singu)),'o-');
%~ xlabel('Lcav');
%~ ylabel(['module de l admittance à f=' num2str(f)]);

[s pos]=max(admittance_singu);


disp(['Avec la formule, Lcav = ' num2str(Lcav) ' m. Par calcul matriciel : ' num2str(L(pos)) ' m'])


endfunction
