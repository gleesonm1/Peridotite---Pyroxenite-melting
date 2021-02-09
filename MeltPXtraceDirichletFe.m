function [Ctfinal,prope,CREEf,CREE,CfH2O,Hflux,Urin,Pfinal,Cmix,CH2OGarnet,Feiso,CFlourine]=MeltPXtraceDirichletFe(Pressure,F_Pyx,F_Per,T,Dsource,Esource,dtop,Ur,Uchange,Dchange,Xd,Ffactor,m,N,dFePdstart,dFePystart)
rho=3300;
g=9.81;

%% calculate mineral-melt and mineral-mineral fractionation factors
Echarge=1.60217662*10^(-19); % charge of an electron
Econst=8.85418782*10^(-12); %electric constant
Fe3Fet_meltPd=0.13;
Fe3Fet_meltPy=0.25;
Fe3Fet_cpx=0.055;
Fe3Fet_opx=0.035;
Fe3Fet_grt=0.05;
nnn=12;

z_O=2;
z_ol=2;
z_spl=2;
z_cpx=Fe3Fet_cpx*3+(1-Fe3Fet_cpx)*2;
z_opx=Fe3Fet_opx*3+(1-Fe3Fet_opx)*2;
z_grt=Fe3Fet_grt*3+(1-Fe3Fet_grt)*2;

% calculate the force constants for the different minerals
% % KfOl=(2*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.38+0.78)*10^(-10))^3));
KfOl=-(6/(4*0.36))*(((2/6)*(2/4)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.38+0.78)*10^(-10))^3));
% % KfCpx=(1/3)*(Fe3Fet_cpx*((3*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.36+0.645)*10^(-10))^3)))+(1-Fe3Fet_cpx)*((2*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.36+0.78)*10^(-10))^3))))+...
% %     (2/3)*(Fe3Fet_cpx*((3*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.38+0.645)*10^(-10))^3)))+(1-Fe3Fet_cpx)*((2*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.38+0.78)*10^(-10))^3))));
KfCpx=-(6/(4*0.36))*((1/3)*((1-Fe3Fet_cpx)*(((2/6)*(2/3)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.36+0.78)*10^(-10))^3))+Fe3Fet_cpx*(((3/6)*(2/3)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.36+0.645)*10^(-10))^3)))+...
    (2/3)*((1-Fe3Fet_cpx)*(((2/6)*(2/4)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.38+0.78)*10^(-10))^3))+Fe3Fet_cpx*(((3/6)*(2/4)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.38+0.645)*10^(-10))^3))));
% % KfSp=(2*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.38+0.78)*10^(-10))^3));
KfSp=-(4/(4*0.36))*(((2/4)*(2/4)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.38+0.63)*10^(-10))^3));
yPd=(351-199)*Fe3Fet_meltPd+199;
yPy=(351-199)*Fe3Fet_meltPy+199;
% % KfOpx=0.5*((1/3)*(2*Fe3Fet_opx*((3*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.36+0.645)*10^(-10))^3)))+(1-2*Fe3Fet_opx)*((2*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.36+0.78)*10^(-10))^3))))+...
% %     (2/3)*(2*Fe3Fet_opx*((3*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.38+0.645)*10^(-10))^3)))+(1-2*Fe3Fet_opx)*((2*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.38+0.78)*10^(-10))^3)))))+...
% %     0.5*(0.5*((3/7)*(((2*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.36+0.78)*10^(-10))^3))))+(4/7)*(((2*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.38+0.78)*10^(-10))^3)))))+...
% %     0.5*((3/7)*(((2*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.36+0.92)*10^(-10))^3))))+(4/7)*(((2*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.38+0.92)*10^(-10))^3))))));
KfOpx=-(0.5*(6/(4*0.36))*((1/3)*(Fe3Fet_opx*(((3/6)*(2/3)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.36+0.645)*10^(-10))^3))+...
    (1-Fe3Fet_opx)*(((2/6)*(2/3)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.36+0.78)*10^(-10))^3)))+(2/3)*(Fe3Fet_opx*(((3/6)*(2/4)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.38+0.645)*10^(-10))^3))+...
    (1-Fe3Fet_opx)*(((2/6)*(2/4)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.38+0.78)*10^(-10))^3))))+0.5*(0.5*(6/(4*0.36))*((3/7)*(Fe3Fet_opx*(((3/6)*(2/3)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.36+0.645)*10^(-10))^3))+...
    (1-Fe3Fet_opx)*(((2/6)*(2/3)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.36+0.78)*10^(-10))^3)))+(4/7)*(Fe3Fet_opx*(((3/6)*(2/4)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.38+0.645)*10^(-10))^3))+...
    (1-Fe3Fet_opx)*(((2/6)*(2/4)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.38+0.78)*10^(-10))^3))))+0.5*(8/(4*0.36))*((3/7)*(Fe3Fet_opx*(((3/8)*(2/3)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.36+0.78)*10^(-10))^3))+...
    (1-Fe3Fet_opx)*(((2/8)*(2/3)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.36+0.92)*10^(-10))^3)))+(4/7)*(Fe3Fet_opx*(((3/8)*(2/4)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.38+0.78)*10^(-10))^3))+...
    (1-Fe3Fet_opx)*(((2/8)*(2/4)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.38+0.92)*10^(-10))^3))))));
% % KfGrt=Fe3Fet_grt*(3*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.38+0.645)*10^(-10))^3))+(1-Fe3Fet_grt)*(2*(-2)*(Echarge^2)*(1-nnn))/(4*pi*Econst*(((1.38+0.92)*10^(-10))^3));
KfGrt=-(Fe3Fet_grt*((6/(4*0.36))*(((3/6)*(2/4)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.38+0.645)*10^(-10))^3)))+(1-Fe3Fet_grt)*((8/(4*0.36))*(((2/8)*(2/4)*((1.602*10^(-19))^2)*(1-12)*(9*10^9))/(((1.38+0.92)*10^(-10))^3))));

%% Load fractions vs depth
Pressure=(round(Pressure.*100))./100;
dh=((9-8.99)*10^9)/(rho*g); % change of z per pressure step in melt px model

AAA=find(F_Pyx>0,1); % find the point at which pyroxenite starts to melt
BBB=find(F_Per>0,1); % find the point at which peridotite starts to melt

k=find(Pressure==6); % find point at which P=6GPa
top=find(Pressure==dtop); % find top of melt column
final=[Pressure F_Pyx F_Per];
final=final(AAA-1:top,:); % only take forward value between the start of pyroxenite melting and the top of the melt column
AAA=find(final(:,2)>0,1); % new location for start of pyroxenite melting
BBB=find(final(:,3)>0,1); % new location for start of peridotite melting

Zo=((final(AAA-1,1)-dtop)*10^9)/(rho*g);  % depth to the start of pyroxenite melting
ZoD=((final(BBB-1,1)-dtop)*10^9)/(rho*g); % depth to the start of peridotite melting
TF=(T(AAA-1:top));

final(:,2)=final(:,2)./100;
final(:,3)=final(:,3)./100;

%% calculate Fe-isotope fractionation factors
Pd_ol=2904.*(yPd-KfOl)./((TF+273.15).^2);%-(2./3).*(-0.4).*(10.^6)./((TF+273.15).^2);%2853.*(yPd-133)./((TF+273.15).^2); %KfOl
Py_ol=2904.*(yPy-KfOl)./((TF+273.15).^2);%-(2./3).*(-0.4).*(10.^6)./((TF+273.15).^2);%2853.*(yPy-133)./((TF+273.15).^2);
ol_cpx=2904.*(KfOl-KfCpx)./((TF+273.15).^2);
ol_opx=2904.*(KfOl-KfOpx)./((TF+273.15).^2);
ol_sp=2904.*(KfOl-KfSp)./((TF+273.15).^2);
ol_grt=2904.*(KfOl-KfGrt)./((TF+273.15).^2);

DFeOl=1; DFeOpx=0.6; DFeCpx=0.3; DFeSp=2; DFeGrt=1.6;
nOl=0.57; nOpx=0.26; nCpx=0.15; nSp=0.02; nPyCpx=0.7; nPyGrt=0.3;
FeOl=10.5; FeOpx=6; FeCpx=3; FeSp=15; FePyCpx=5.9; FePyGrt=14;
pOl=-0.2; pOpx=0.61; pCpx=0.58; pSp=0.01; pPyCpx=0.9; pPyGrt=0.1;
%% define volumes
% Zo=((6-dtop)*10^9)/(rho*g);
dc=((Dchange-dtop)*10^9)/(rho*g); 

for i=1:length(final(:,1))
    if final(i,1)>Dchange
        v(i)=(2*dh)/(dc+2*(Zo-dc));
    else
        n=i-find(final(:,1)==Dchange)+1;
        v(i)=(2*dh*(dc-(n-0.5)*dh))/(dc*(dc+2*(Zo-dc)));
    end
end
    


%% Define partition coefficients
% olivine
DBaOl=0.000005; DNbOl=0.0005; DPbOl=0.003; DSrOl=0.00004; DLaOl=0.0005;
DCeOl=0.0005; DPrOl=0.0008; DNdOl=0.00042; DSmOl=0.0011; DEuOl=0.0016;
DGdOl=0.0011; DTbOl=0.0015; DDyOl=0.0027; DHoOl=0.0016; DErOl=0.013; 
DYbOl=0.02; DLuOl=0.02; DH2OOl=0.00125; DFOl=0.005;
DOl=[DPbOl DSrOl DLaOl DCeOl DPrOl DNdOl DSmOl DEuOl DGdOl DTbOl DDyOl DHoOl DErOl DYbOl DLuOl DH2OOl DFOl];

% orthopyroxene
DBaOp=0.000006; DNbOp=0.004; DPbOp=0.009; DSrOp=0.0007; DLaOp=0.0031;
DCeOp=0.004; DPrOp=0.0048; DNdOp=0.012; DSmOp=0.02; DEuOp=0.013;
DGdOp=0.0065; DTbOp=0.019; DDyOp=0.011; DHoOp=0.026; DErOp=0.045; 
DYbOp=0.08; DLuOp=0.12; DH2OOp=0.0145; DFOp=0.05;
DOp=[DPbOp DSrOp DLaOp DCeOp DPrOp DNdOp DSmOp DEuOp DGdOp DTbOp DDyOp DHoOp DErOp DYbOp DLuOp DH2OOp DFOp];

% Garnet
DBaG=0.00007; DNbG=0.015; DPbG=0.005; DSrG=0.0007; DLaG=0.001;
DCeG=0.005; DPrG=0.014; DNdG=0.052; DSmG=0.25; DEuG=0.496;
DGdG=0.848; DTbG=1.477; DDyG=2.2; DHoG=3.315; DErG=4.4; 
DYbG=6.6; DLuG=7.1; DH2OG=0.00316; DFG=0.0123;
DG=[DPbG DSrG DLaG DCeG DPrG DNdG DSmG DEuG DGdG DTbG DDyG DHoG DErG DYbG DLuG DH2OG DFG];

% clinopyroxene
DBaCp=0.0004; DNbCp=0.015; DPbCp=0.012; DSrCp=0.091; DLaCp=0.049;
DCeCp=0.08; DPrCp=0.126; DNdCp=0.178; DSmCp=0.293; DEuCp=0.335;
DGdCp=0.35; DTbCp=0.403; DDyCp=0.4; DHoCp=0.427; DErCp=0.42; 
DYbCp=0.4; DLuCp=0.376; DH2OCp=0.0139; DFCp=0.13;
DCp=[DPbCp DSrCp DLaCp DCeCp DPrCp DNdCp DSmCp DEuCp DGdCp DTbCp DDyCp DHoCp DErCp DYbCp DLuCp DH2OCp DFCp];

% spinel
DBaS=0.000; DNbS=0.0; DPbS=0.0; DSrS=0.0; DLaS=0.0006;
DCeS=0.0006; DPrS=0.0006; DNdS=0.0006; DSmS=0.0006; DEuS=0.0006;
DGdS=0.0006; DTbS=0.0006; DDyS=0.0015; DHoS=0.0023; DErS=0.003; 
DYbS=0.0045; DLuS=0.0052; DH2OS=0.0006; DFS=0.0006;
DS=[DPbS DSrS DLaS DCeS DPrS DNdS DSmS DEuS DGdS DTbS DDyS DHoS DErS DYbS DLuS DH2OS DFS];

%% define compositions

if Dsource==1
    % composition of WH05 DDM with ~2% melt extraction in the garnet
    % stability field
    Cdin=[0.011 5.33478 0.092921 0.343447 0.079655 0.471929 0.21581 0.089761 0.343601 0.068199 0.495353 0.113472 0.344595 0.364105 0.057964 65 9];
end

if Dsource==2
    % composition of WH05 DDM with ~2% melt extraction in the garnet
    % stability field combined with the Donelley et al. 2004 estimate for
    % the enriched mantle in a 0.95:0.05 ratio 0.9:0.1
    % Cdin=[0.01705 5.993541 0.125125 0.408775 0.088172 0.509833 0.22462 0.092723 0.353221 0.069639 0.503885 0.115198 0.349215 0.364795 0.058466 67.5 11.05];
    Cdin=[0.0201 6.65302 0.157329 0.474102 0.096689 0.547736 0.233429 0.095685 0.362841 0.071079 0.512418 0.116925 0.353836 0.371795 0.058968 95 13.1];
end

if Dsource==3
    % Depleted DMM composition from Workman and Hart 2005
    Cdin=[0.014 6.092 0.134 0.421 0.087 0.483 0.21 0.086 0.324 0.064 0.471 0.108 0.329 0.348 0.056 100 12];
end

if Dsource==4
    % Depleted DMM composition from Workman and Hart 2005
    Cdin=[0.0201 7.3338 0.1943 0.5439 0.1033 0.5577 0.2282 0.0923 0.3452 0.0673 0.4905 0.112 0.3398 0.3573 0.0572 141 15];
end

if Esource==1
    % Primitve mantle composition from Sun and Mcdonough 1989
    Cpin=[0.071 21.1 0.687 1.775 0.276 1.354 0.444 0.168 0.596 0.108 0.737 0.164 0.48 0.493 0.074 280 36];
end

if Esource==2
    % 50:50 WHDDM with average subducted MORB from Porter and White 2009
    Cpin=[0.3355 66.05 2.7285 6.9375 1.078 5.267 1.667 0.614 2.148 0.364 2.5935 0.507 1.59 1.5265 0.23 950 115];
end

if Esource==3
    % 50:50 WHDDM with average subducted MORB from Porter and White 2009
    Cpin=[0.16 52.101 1.415 4.291 0.6677 4.429 1.26 0.571 2.012 0.382 2.436 0.629 1.491 1.874 0.214 677 80];
end

if Esource==4
    % 50:50 WHDDM with average subducted MORB from Stracke et al. 2003
    Cpin=[0.052 43.546 0.907 3.1555 0.7435 3.9665 1.45 0.563 2.177 0.342 2.7405 0.479 1.7295 1.669 0.253 350 70];
end

%% melting equation
% spinel peridotite melting equations
% Prop ol 0.578 opx 0.27 cpx 0.119 spl 0.033
P=DCp.*0.5+DOl.*0.1+DOp.*0.27+DS.*0.13;

cpx=0.119-final(:,3).*0.5; ol=0.578-final(:,3).*0.1; opx=0.27-final(:,3).*0.27; spl=0.033-final(:,3).*0.13;
cpx1=cpx./(cpx+ol+opx+spl); ol1=ol./(cpx+ol+opx+spl); opx1=opx./(cpx+ol+opx+spl); spl1=spl./(cpx+ol+opx+spl);

D=cpx1.*DCp+ol1.*DOl+opx1.*DOp+spl1.*DS;

for i=1:length(final(:,3))
    if i==1
        F(i)=final(i,3);
        Fmass(i)=final(i,3);
    else
        F(i)=(final(i,3)-final(i-1,3))/(1-final(i-1,3));
        Fmass(i)=final(i,3)-final(i-1,3);
    end
end

% F=F';


for i=1:length(F)
    if final(i,1)<Uchange
        Urin(i)=1;
    else
        Urin(i)=1+Ur*(2*((final(i,1)-Uchange)/(final(1,1)-Uchange))-((final(i,1)-Uchange)/(final(1,1)-Uchange))^2);
    end
end

% C=zeros(length(final(:,3)),1);
for i=1:length(final(:,3))
    if final(i,3)>0
    if i==BBB
        C(i,:)=Cdin./(D(i,:)+F(i).*(1-P));
        C(isnan(C))=0;
        Res=Cdin.*((D(i,:)-P.*F(i))/(1-F(i))).*(1./(D(i,:)+F(i).*(1-P)));
        
        % calculate the Fe-isotope composition
        nOl=nOl-F(i)*pOl; nOpx=nOpx-F(i)*pOpx; nCpx=nCpx-F(i)*pCpx; nSp=nSp-F(i)*pSp;
        FeMelt(i)=8.3./(F(i)+(nOl*DFeOl+nCpx*DFeCpx+nOpx*DFeOpx+nSp*DFeSp));
        FeOl=FeMelt(i)*DFeOl; FeCpx=FeMelt(i)*DFeCpx; FeOpx=FeMelt(i)*DFeOpx; FeSp=FeMelt(i)*DFeSp; 
        d56FeMelt(i)=(dFePdstart*8.3+Pd_ol(i)*(nOl*FeOl+nCpx*FeCpx+nOpx*FeOpx+nSp*FeSp)+(nCpx*FeCpx*(ol_cpx(i))+nOpx*FeOpx*(ol_opx(i))+nSp*FeSp*(ol_sp(i))))/(F(i)*FeMelt(i)+nOl*FeOl+nCpx*FeCpx+nOpx*FeOpx+nSp*FeSp);
        nOl=nOl/(nOl+nCpx+nOpx+nSp); nCpx=nCpx/(nOl+nCpx+nOpx+nSp); nOpx=nOpx/(nOl+nCpx+nOpx+nSp); nSp=nSp/(nOl+nCpx+nOpx+nSp);
        FeRes(i)=nOl*FeOl+nCpx*FeCpx+nOpx*FeOpx+nSp*FeSp;
        d56FeRes(i)=(dFePdstart*8.3-(F(i)*d56FeMelt(i)*FeMelt(i)))/((1-F(i))*FeRes(i));        
    else
        C(i,:)=Res./(D(i,:)+F(i).*(1-P));
        C(isnan(C))=0;
        Res=Res.*((D(i,:)-P.*F(i))/(1-F(i))).*(1./(D(i,:)+F(i).*(1-P)));
        
        % calculate the Fe-isotope composition
        nOl=nOl-F(i)*pOl; nOpx=nOpx-F(i)*pOpx; nCpx=nCpx-F(i)*pCpx; nSp=nSp-F(i)*pSp;
        FeMelt(i)=FeRes(i-1)./(F(i)+(nOl*DFeOl+nCpx*DFeCpx+nOpx*DFeOpx+nSp*DFeSp));
        FeOl=FeMelt(i)*DFeOl; FeCpx=FeMelt(i)*DFeCpx; FeOpx=FeMelt(i)*DFeOpx; FeSp=FeMelt(i)*DFeSp; 
        d56FeMelt(i)=(d56FeRes(i-1)*FeRes(i-1)+Pd_ol(i)*(nOl*FeOl+nCpx*FeCpx+nOpx*FeOpx+nSp*FeSp)+(nCpx*FeCpx*(ol_cpx(i))+nOpx*FeOpx*(ol_opx(i))+nSp*FeSp*(ol_sp(i))))/(F(i)*FeMelt(i)+nOl*FeOl+nCpx*FeCpx+nOpx*FeOpx+nSp*FeSp);
        nOl=nOl/(nOl+nCpx+nOpx+nSp); nCpx=nCpx/(nOl+nCpx+nOpx+nSp); nOpx=nOpx/(nOl+nCpx+nOpx+nSp); nSp=nSp/(nOl+nCpx+nOpx+nSp);
        FeRes(i)=nOl*FeOl+nCpx*FeCpx+nOpx*FeOpx+nSp*FeSp;
        d56FeRes(i)=(d56FeRes(i-1)*FeRes(i-1)-(F(i)*d56FeMelt(i)*FeMelt(i)))/((1-F(i))*FeRes(i)); 
    end
    end
end
C(isnan(C))=0;

for i=1:length(Fmass)
    Frac(i)=Fmass(i).*Urin(i).*v(i);
end
Frac=Frac./sum(Frac);

for i=1:length(Fmass)
    Cf(i,:)=C(i,:).*Frac(i);
    d56FePd(i)=d56FeMelt(i)*FeMelt(i)*Frac(i);
    FePd(i)=FeMelt(i)*Frac(i);
end
FePdF=sum(FePd);
d56FePdF=sum(d56FePd)/FePdF;
for i=1:length(D(1,:))
    Cfinal(i)=sum(Cf(:,i));
end

%% enriched melt

cpx=0.7-final(:,2).*0.9; grt=0.3-final(:,2).*0.1;
for i=1:length(cpx)
    if cpx(i)>0
        cpx1(i)=cpx(i)./(cpx(i)+grt(i)); grt1(i)=grt(i)./(cpx(i)+grt(i));
        PE(i,:)=DG.*0.1+DCp.*0.9;
    else
        cpx1(i)=0; grt1(i)=1;
        PE(i,:)=DG;
    end
end

DE=cpx1.*DCp+grt1'.*DG;

for i=1:length(final(:,2))
    if i==1
        Fe(i)=final(i,2);
        Femass(i)=final(i,2);
    else
        Fe(i)=(final(i,2)-final(i-1,2))/(1-final(i-1,2));
        Femass(i)=final(i,2)-final(i-1,2);
    end
end

Fe=Fe';

for i=1:length(final(:,3))
    if i==1
        CE(i,:)=Cpin./(DE(i,:)+Fe(i).*(1-PE(i,:)));
        CE(isnan(CE))=0;
        ResE=Cpin.*((DE(i,:)-PE(i,:).*Fe(i))/(1-Fe(i))).*(1./(DE(i,:)+Fe(i).*(1-PE(i,:))));
        
        % calculate the Fe-isotope composition
        nPyCpx=nPyCpx-Fe(i)*pPyCpx; nPyGrt=nPyGrt-Fe(i)*pPyGrt;
        FeMeltPy(i)=8.3/(Fe(i)+nPyGrt*DFeGrt+nPyCpx*DFeCpx);
        FePyCpx=FeMeltPy(i)*DFeCpx; FePyGrt=FeMeltPy(i)*DFeGrt;
        d56FeMeltPy(i)=(dFePystart*8.3+Py_ol(i)*(nPyGrt*FePyGrt+nPyCpx*FePyCpx)+(nPyGrt*FePyGrt*(ol_grt(i))+nPyCpx*FePyCpx*(ol_cpx(i))))/(Fe(i)*FeMeltPy(i)+nPyGrt*FePyGrt+nPyCpx*FePyCpx);
        nPyGrt=nPyGrt/(nPyGrt+nPyCpx); nPyCpx=nPyCpx/(nPyGrt+nPyCpx);
        FeResPy(i)=nPyCpx*FePyCpx+nPyGrt*FePyGrt;
        d56FeResPy(i)=(dFePystart*8.3-(Fe(i)*d56FeMeltPy(i)*FeMeltPy(i)))/((1-Fe(i))*FeResPy(i));

    else
        CE(i,:)=ResE./(DE(i,:)+Fe(i).*(1-PE(i,:)));
        CE(isnan(CE))=0;
        ResE=ResE.*((DE(i,:)-PE(i,:).*Fe(i))/(1-Fe(i))).*(1./(DE(i,:)+Fe(i).*(1-PE(i,:))));
        
        % calculate the Fe-isotope composition
        nPyCpx=nPyCpx-Fe(i)*pPyCpx; nPyGrt=nPyGrt-Fe(i)*pPyGrt;
        FeMeltPy(i)=FeResPy(i-1)/(Fe(i)+nPyGrt*DFeGrt+nPyCpx*DFeCpx);
        FePyCpx=FeMeltPy(i)*DFeCpx; FePyGrt=FeMeltPy(i)*DFeGrt;
        d56FeMeltPy(i)=(d56FeResPy(i-1)*FeResPy(i-1)+Py_ol(i)*(nPyGrt*FePyGrt+nPyCpx*FePyCpx)+(nPyGrt*FePyGrt*(ol_grt(i))+nPyCpx*FePyCpx*(ol_cpx(i))))/(Fe(i)*FeMeltPy(i)+nPyGrt*FePyGrt+nPyCpx*FePyCpx);
        nPyGrt=nPyGrt/(nPyGrt+nPyCpx); nPyCpx=nPyCpx/(nPyGrt+nPyCpx);
        FeResPy(i)=nPyCpx*FePyCpx+nPyGrt*FePyGrt;
        d56FeResPy(i)=(d56FeResPy(i-1)*FeResPy(i-1)-(Fe(i)*d56FeMeltPy(i)*FeMeltPy(i)))/((1-Fe(i))*FeResPy(i));
    end
end
CE(isnan(CE))=0;

for i=1:length(Femass)
    FracE(i)=Femass(i).*Urin(i).*v(i);
end
FracE=FracE./sum(FracE);

for i=1:length(Femass)
    Cfe(i,:)=CE(i,:).*FracE(i);
    d56FePy(i)=d56FeMeltPy(i)*FeMeltPy(i)*FracE(i);
    FePy(i)=FeMeltPy(i)*FracE(i);
end
FePyF=sum(FePy);
d56FePyF=sum(d56FePy)/FePyF;
for i=1:length(DE(1,:))
    Cfinale(i)=sum(Cfe(:,i));
end


%% calculating crustal thickness
% need to calculate heights and B factors already have dc, but now need to
% redesign Zo and calculate duc
duc=((Uchange-dtop)*10^9)/(rho*g);
Zd=find(final(:,3),1); Zod=((final(Zd-1,1)-dtop)*10^9)/(rho*g);
Ze=find(final(:,2),1); Zoe=((final(Ze-1,1)-dtop)*10^9)/(rho*g);

for i=1:length(final(:,1))
    Bpd(i)=final(i,3)./(Zod-((final(i,1)-dtop)*10^9)/(rho*g));
    Bpx(i)=final(i,2)./(Zoe-((final(i,1)-dtop)*10^9)/(rho*g));
end
Bpd(isnan(Bpd))=0;
Bpx(isnan(Bpx))=0;
UrCT=flip(Urin);
Bpd=flip(Bpd);
Bpx=flip(Bpx);
depth=flip(final(:,1));

for i=1:length(UrCT)
    if i==1
    CTd(i)=UrCT(i)*((Bpd(i)*(Zod^2))/2);
    CTe(i)=UrCT(i)*((Bpx(i)*(Zoe^2))/2);
    else
        CTd(i)=(UrCT(i)-UrCT(i-1))*((Bpd(i)*((Zod^2)-((((depth(i)-dtop)*10^9)/(rho*g))^2)))/2);
        CTe(i)=(UrCT(i)-UrCT(i-1))*((Bpx(i)*((Zoe^2)-((((depth(i)-dtop)*10^9)/(rho*g))^2)))/2);
    end
end
CTD=sum(CTd);
CTE=sum(CTe);
TT=find(depth==Dchange);




% calculated hypothetical thickness below Dchange
for i=TT:length(UrCT)
    if i==TT
    CTd2(i)=UrCT(i)*((Bpd(i)*((Zod^2)-((((depth(i)-dtop)*10^9)/(rho*g))^2)))/2);
    CTe2(i)=UrCT(i)*((Bpx(i)*((Zoe^2)-((((depth(i)-dtop)*10^9)/(rho*g))^2)))/2);
    else
        CTd2(i)=(UrCT(i)-UrCT(i-1))*((Bpd(i)*((Zod^2)-((((depth(i)-dtop)*10^9)/(rho*g))^2)))/2);
        CTe2(i)=(UrCT(i)-UrCT(i-1))*((Bpx(i)*((Zoe^2)-((((depth(i)-dtop)*10^9)/(rho*g))^2)))/2);
    end
end
CTD2=sum(CTd2);
CTE2=sum(CTe2);

CTDF=CTD-CTD2;
CTEF=CTE-CTE2;

% calculate proportion of melt from above Dchange

Fd=flip(Frac);
Fe=flip(FracE);
PropAd=sum(Fd(1:TT));
PropAe=sum(Fe(1:TT));

CTde=(3300/2900)*CTDF/(PropAd);
CTee=(3400/2900)*CTEF/PropAe;

Ctfinal=(Xd*CTde+(1-Xd)*CTee);

%% Calculate final melt composition

% Pdd=CTde/Ctfinal;
Pdd=sum(Urin.*v.*Fmass.*Xd)/(sum(Urin.*v.*Fmass.*Xd)+sum(Urin.*v.*Femass.*(1-Xd)));
prope=(1-Pdd)*(1-m)+m;

Feiso=(Pdd*FePdF*d56FePdF+(1-Pdd)*FePyF*d56FePyF)/(Pdd*FePdF+(1-Pdd)*FePyF);
Cff=Pdd.*Cfinal+(1-Pdd).*Cfinale;

PPb=1; PSr=2; PLa=3; PCe=4; PPr=5; PNd=6; PSm=7; PEu=8; PGd=9; PTb=10; PDy=11; PHo=12; PEr=13; PYb=14; PLu=15; PH=16; PF=17;

% CfH2O=Pdd.*Cfinal(PH)+(1-Pdd).*Cfinale(PH);


% Srloc=find(strcmp(raw_P(1,:),'Sr')); Pbloc=find(strcmp(raw_P(1,:),'Pb'));
% CSrd=Cfinal(Srloc); CSre=Cfinale(Srloc); CSrf=Cff(Srloc);
% CPbd=Cfinal(Pbloc); CPbe=Cfinale(Pbloc); CPbf=Cff(Pbloc);

% Pb208=(CPbd*Pdd*EMiso(2,1)+CPbe*(1-Pdd)*EMiso(4,1))/CPbf;
% Pb207=(CPbd*Pdd*EMiso(2,2)+CPbe*(1-Pdd)*EMiso(4,2))/CPbf;
% Pb206=(CPbd*Pdd*EMiso(2,3)+CPbe*(1-Pdd)*EMiso(4,3))/CPbf;
% Sr87=(CSrd*Pdd*EMiso(2,4)+CSre*(1-Pdd)*EMiso(4,4))/CSrf;

% norm=[La Cerium Pr Nd Sm Eu Gd Tb Dy Ho Er Yb Lu];
% if choice==2
% Cgarnet=[20.53506746 8714 48.60205853 5.722630807 26.2043574 40 5.382610145 1.58807228 4.646289438 0.626585627 3.453266796 0.569395012 1.281145681 0.910924834 0.13019392];
% elseif choice==1
%     Cgarnet=[54.39432575 24618.82326 109.8059859 12.4309263 40.97852032 40 6.443361413 1.696964094 4.344650751 0.504894183 2.532301172 0.416275416 0.881438505 0.622448516 0.086883653];
% end
Cgarnet=(Cpin./0.10).*(1-(1-(PE(1,:).*0.1)./DE(1,:)).^(1./PE(1,:)));

CffG=m.*Cgarnet+(1-m).*Cff;    
CREE=[CffG(PLa) CffG(PCe) CffG(PPr) CffG(PNd) CffG(PSm) CffG(PEu) CffG(PGd) CffG(PTb) CffG(PDy) CffG(PHo) CffG(PEr) CffG(PYb) CffG(PLu)];

CfH2O=CffG(PH);
CFlourine=CffG(PF);
Ctfinal=Ctfinal./(1-m);
Hflux=(CffG(PH)*10^(-6))*0.057*Ctfinal*2900;

CH2OGarnet=Cgarnet(PH);
CREEf=CREE./Ffactor;
Pfinal=final(:,1);

%% working out Dirichlet mixing function
for i=1:length(C(:,1))
Pmd(i)=(Urin(i).*v(i).*Fmass(i).*Xd)/(sum(Urin.*v.*Fmass.*Xd)+sum(Urin.*v.*Femass.*(1-Xd)));
end
for i=1:length(CE(:,1))
    Pme(i)=(Urin(i).*v(i).*Femass(i).*(1-Xd))/(sum(Urin.*v.*Fmass.*Xd)+sum(Urin.*v.*Femass.*(1-Xd)));
end
% for ii=1:length(Frac)
%     if final(ii,3)==0
%         Frac(ii)=0;
%     end
% end
% Frac=Frac./sum(Frac);
% Pmd=Frac.*Pdd;
% Pme=FracE.*(1-Pdd);
CDiri=[C;CE;Cgarnet];
Pmass=[Pmd Pme]';
Pmass=[Pmass*(1-m);m];
for ii=1:length(Pmass)
    if Pmass(ii)<0
        Pmass(ii)=0;
    end
end

for iii=1:500

    alpha=(N-1).*Pmass;
    r=dirichletRnd(alpha',1);
    for ii=1:length(CDiri(1,:))
    Cmix(iii,ii)=sum(r'.*CDiri(:,ii));
    end
end
end