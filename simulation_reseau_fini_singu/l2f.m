function l2f(Lcav) %Lcol est la longueur de la cavité

Lcol=0.02; %longueur du col
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
Lcol = Lcol + L1 + L2;
c=343;

%Par calcul matriciel
f=0:1:1000;
w=2*pi*f;
N=length(f);

for x=1:1:N
	matrice_resonateur = resonateur(w(x),Lcav,Lcol,Dcav,Dcol,1.2,c);		
	admittance_singu(x) = matrice_resonateur(2,1);
end

[s pos]=max(admittance_singu);
f1=f(pos);


%Avec la formule
f2=c/(2*pi) *sqrt(Scol/(Scav*Lcav*Lcol));

disp(['Avec  la formule : f = ' num2str(f2) ' Hz. Par le calcul matriciel : ' num2str(f1) ])

endfunction
