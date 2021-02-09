function [Pressure,F_Per,F_Pyx,T,T_solidadiabat,TsolPer,TsolPyx]=MeltPX(MantleTp,FractionPyx,FractionCpx,inSiO2,inTiO2,inAl2O3,inCr2O3,inFeO,inMnO,inMgO,inCaO,inNa2O,inK2O)
    
%% set initial parameters
% user input parameters
Tp=MantleTp; % oC mantle potential temperature
FracPyx=FractionPyx; % Fraction of pyroxenite in the mantle source
% FracPer=(1-FracPyx); % Fraction of peridotite in the mantle source
FracCpx=FractionCpx; % Mass fraction of cpx in the subsolidus peridotite

% composition of pyroxenite (in mol%)
SiO2=inSiO2; %47.5395553422756;
TiO2=inTiO2; %0.3832220670177;
Al2O3=inAl2O3; %7.14067512672454;
Cr2O3=inCr2O3; %0.0464700082784867;
FeO=inFeO; %7.38957897155817;
MnO=inMnO; %0.165943726786654;
MgO=inMgO; %24.3001432174423;
CaO=inCaO; %11.4300141859886;
Na2O=inNa2O; %1.56690635075243;
K2O=inK2O; %0.0374910031755343;
MgNo=MgO/(MgO+FeO);

%% pre defined parameters
% pyroxenite
CpPy=1140; % hear capacity
alphaPy=0.00003; % thermal expansivity
rhoPy=3400; % density of solid
SPy=380; %Entropy of fusion

% peridotite
CpPd=1190; % hear capacity
alphaPd=0.00003; % thermal expansivity
rhoPd=3300; % density of solid
SPd=395; %Entropy of fusion

% Melt-PX parameters
a=1.4; b=-7.08; c=-12.5; d=150; e=25.5; f=1105; % a, b, c, d, e, f
aa=-14; bb=177; cc=219; dd=-17.2; ee=10.4; ff=1000; % a', b', c', d', e', f'
A=-2.63; B=-1.58; C=26.9; D=-7.5; E=4.23; AA=2.86; BB=18; % A', B', C', D', E', A'', B''

% density of liquid
rhoL=2900;

%% calculate P and T of solid adiabat
P=linspace(6,0.20,581);
P=P';

% solid adiabat
T_SA=Tp+((1-FracPyx).*((alphaPd.*(Tp+273.15))./(CpPd.*rhoPd).*10^9).*P+FracPyx.*((alphaPy.*(Tp+273.15))./(CpPy.*rhoPy).*10^9).*P);

% Paramaterisation for pyroxenite
% T 5% pyroxenite
T5Pyx=(a.*(Na2O+K2O)+b).*(P(1:end-1).^2)+(c.*(Na2O+K2O)+d.*MgNo+e).*P(1:end-1)+f;

% T cpx-out pyroxenite
Tcpx_outPyx=aa.*(P(1:end-1).^2)+bb.*P(1:end-1)+cc.*MgNo+dd.*Al2O3+ee.*CaO+ff;

%T' pyroxenite
TprimPyx=(-(AA.*CaO+BB)+((AA.*CaO+BB).^2-4.*5.*((A+B.*MgNo).*(P(1:end-1).^2)+(C+D.*MgNo).*P(1:end-1)+E)).^0.5)./(2.*((A+B.*MgNo).*(P(1:end-1).^2)+(C+D.*MgNo).*P(1:end-1)+E));

%T_sol pyroxenite
TsolPyx=TprimPyx.*(Tcpx_outPyx-T5Pyx)+T5Pyx;

% Parameterisation for peridotite
%T_sol peridotite
TsolPer=1085.7+132.9.*P(1:end-1)-5.1.*(P(1:end-1).^2);

% T_lherz^liq peridotite
Tliqlherz=1475+80.*P(1:end-1)-3.2.*P(1:end-1).^2;

% Fcpx-out
Fcpx_out=FracCpx./(0.5+0.08.*P(1:end-1));

% Tcpx-out peridotite
Tcpx_outPer=(Fcpx_out.^(1./1.5)).*(Tliqlherz-TsolPer)+TsolPer;

% Tliquidus - peridotite
TliqPer=1780+45.*P(1:end-1)-2.*P(1:end-1).^2;

% T'(Pyx)_solid adiabat
TprimPyx_SA=zeros(length(P)-1,1);
for i=1:length(P)-1
    if Tcpx_outPyx(i)<T5Pyx(i)
        TprimPyx_SA(i)=-1;
    else
        TprimPyx_SA(i)=(T_SA(i)-T5Pyx(i))/(Tcpx_outPyx(i)-T5Pyx(i));
    end
end

% F(Pyx) solid adiabat
FPyx_SA=zeros(length(P)-1,1);
for i=1:length(P)-1
    if TprimPyx_SA(i)<TprimPyx(i)
        FPyx_SA(i)=0;
    else
        FPyx_SA(i)=((A+B*MgNo)*(P(i)^2)+(C+D*MgNo)*P(i)+E)*(TprimPyx_SA(i)^2)+(AA*CaO+BB)*TprimPyx_SA(i)+5;
    end
end


% T' (per) solid adiabat
TprimPer_SA=(T_SA(1:end-1)-TsolPer)./(Tliqlherz-TsolPer);

% F cpx solid adiabat peridotite
FcpxPer_SA=zeros(length(P)-1,1);
for i=1:length(P)-1
    if TprimPer_SA(i)<0
        FcpxPer_SA(i)=0;
    else
        FcpxPer_SA(i)=(TprimPer_SA(i)^1.5)*100;
    end
end

% Fopx Per solid adiabat
FopxPer_SA=zeros(length(P)-1,1);
for i=1:length(P)-1
    if T_SA(i)<Tcpx_outPer(i)
        FopxPer_SA(i)=0;
    else
        FopxPer_SA(i)=((Fcpx_out(i)+(1-Fcpx_out(i))*((T_SA(i)-Tcpx_outPer(i))/(TliqPer(i)-Tcpx_outPer(i))))^1.5)*100;
    end
end

% F(per)Solid adiabat
FPer_SA=zeros(length(P)-1,1);
for i=1:length(P)-1
    if FopxPer_SA(i)>0
        FPer_SA(i)=FopxPer_SA(i);
    else
        FPer_SA(i)=FcpxPer_SA(i);
    end
end

FPyx_SA(end+1)=FPyx_SA(end);

%% now need to calculate all the final properties essentially in one go

% dTs/dP)s pyroxenite
dTsdP_s_Pyx=zeros(length(P)-1,1);
for i=2:length(P)-1
    dTsdP_s_Pyx(i)=(TsolPyx(i)-TsolPyx(i-1))/(P(i)-P(i-1));
end
dTsdP_s_Pyx(1)=dTsdP_s_Pyx(2);

% dTs/dP)s peridotite
dTsdP_s_Per=zeros(length(P)-1,1);
for i=2:length(P)-1
    dTsdP_s_Per(i)=(TsolPer(i)-TsolPer(i-1))/(P(i)-P(i-1));
end
dTsdP_s_Per(1)=dTsdP_s_Per(2);

T=zeros(length(P)-1,1);
TprimPyxT=zeros(length(P)-1,1);
FTprimPyx=zeros(length(P)-1,1);
Tpp1_Pyx=zeros(length(P)-1,1);
FTpp1_Pyx=zeros(length(P)-1,1);
dFdT_P_Pyx=zeros(length(P)-1,1);
F_Pyx=zeros(length(P)-1,1);
moddFdT_P_Pyx=zeros(length(P)-1,1);
dTdP_F_Pyx=zeros(length(P)-1,1);
dTdP_F_PyxT=zeros(length(P)-1,1);
dFdP1_Pyx=zeros(length(P)-1,1);
dFdP2_Pyx=zeros(length(P)-1,1);
dFdP3_Pyx=zeros(length(P)-1,1);
dFdP4_Pyx=zeros(length(P)-1,1);
dFdP_S_Pyx=zeros(length(P)-1,1);
TpM_Pyx=zeros(length(P)-2,1);
TM_Pyx=zeros(length(P)-2,1);

TprimPerT=zeros(length(P)-1,1);
F0_cpxT=zeros(length(P)-1,1);
Tpp1_Per=zeros(length(P)-1,1);
F_cpxT=zeros(length(P)-1,1);
FOT=zeros(length(P)-1,1);
FTpp_Per=zeros(length(P)-1,1);
dFdT_P_Per=zeros(length(P)-1,1);
F_Per=zeros(length(P)-1,1);
dTdP_F_Per=zeros(length(P)-1,1);
dTdP_F_PerT=zeros(length(P)-1,1);
dFdP1_Per=zeros(length(P)-1,1);
dFdP2_Per=zeros(length(P)-1,1);
dFdP3_Per=zeros(length(P)-1,1);
dFdP_S_Per=zeros(length(P)-1,1);
TpM_Per=zeros(length(P)-2,1);
TM_Per=zeros(length(P)-2,1);

DTDP=zeros(length(P)-1,1);
dTdP_S_final=zeros(length(P)-1,1);

for i=1:length(P)-1
    if i==1
        T(i)=T_SA(i); % temperature followed by the solid mantle
        if T5Pyx(i)>Tcpx_outPyx(i) % calculating T' for the pyroxenite column
            TprimPyxT(i)=-1;
        else
            TprimPyxT(i)=(T(i)-T5Pyx(i))/(Tcpx_outPyx(i)-T5Pyx(i));
        end
        
        % calculating T' for the peridotite column
        TprimPerT(i)=(T(i)-TsolPer(i))/(Tliqlherz(i)-TsolPer(i));
        
        if TprimPyxT(i)<TprimPyx(i) % F(T') for pyroxenite
            FTprimPyx(i)=0;
        else
            if (((A+B*MgNo)*(P(i)^2)+(C+D*MgNo)*P(i)+E)*(TprimPyxT(i)^2)+(AA*CaO+BB)*TprimPyxT(i)+5)<0
                FTprimPyx(i)=0;
            else
                FTprimPyx(i)=(((A+B*MgNo)*(P(i)^2)+(C+D*MgNo)*P(i)+E)*(TprimPyxT(i)^2)+(AA*CaO+BB)*TprimPyxT(i)+5);
            end
        end
        
        if TprimPerT(i)<0 % F0-CPX (T') for peridotite
            F0_cpxT(i)=0;
        else
            F0_cpxT(i)=(TprimPerT(i)^(1.5))*100;
        end
        
        if T5Pyx(i)>Tcpx_outPyx(i) % T''(T+1) for pyroxenite
            Tpp1_Pyx(i)=-1;
        else
            Tpp1_Pyx(i)=TprimPyxT(i)+(1/(Tcpx_outPyx(i)-T5Pyx(i)));
        end
        
        % T''(T+1) for peridotite
        Tpp1_Per(i)=(T(i)+1-TsolPer(i))/(Tliqlherz(i)-TsolPer(i));
        
        % F-cpx(T'') peridotite
        if Tpp1_Per(i)<0
            F_cpxT(i)=0;
        else
            F_cpxT(i)=(Tpp1_Per(i)^(1.5))*100;
        end
        
        % F0(T') peridotite
        if F0_cpxT(i)<Fcpx_out(i)*100
            FOT(i)=F0_cpxT(i);
        else
            FOT(i)=(Fcpx_out(i)+(1-Fcpx_out(i))*((T(i)-Tcpx_outPer(i))/(TliqPer(i)-Tcpx_outPer(i)))^1.5)*100;
        end
        
        % F(T'') pyroxenite
        if Tpp1_Pyx(i)<TprimPyx(i)
            FTpp1_Pyx(i)=0;
        else
            if (((A+B*MgNo)*(P(i)^2)+(C+D*MgNo)*P(i)+E)*(Tpp1_Pyx(i)^2)+(AA*CaO+BB)*Tpp1_Pyx(i)+5)<0
                FTpp1_Pyx(i)=0;
            else
                FTpp1_Pyx(i)=(((A+B*MgNo)*(P(i)^2)+(C+D*MgNo)*P(i)+E)*(Tpp1_Pyx(i)^2)+(AA*CaO+BB)*Tpp1_Pyx(i)+5);
            end
        end
        
        % F(T'') peridotite
        if F_cpxT(i)<Fcpx_out(i)*100
            FTpp_Per(i)=F_cpxT(i);
        else
            FTpp_Per(i)=(Fcpx_out(i)+(1-Fcpx_out(i))*((T(i)-Tcpx_outPer(i))/(TliqPer(i)-Tcpx_outPer(i)))^1.5);
        end
        
        % Melt fraction of pyroxenitic lithology
        F_Pyx(i)=FPyx_SA(i);
        
        % dF/dT_P for pyroxenite (cpx exhaustion matters)
        if T(i)>Tcpx_outPyx(i)
            dFdT_P_Pyx(i)=0.3;
        else
            dFdT_P_Pyx(i)=FTpp1_Pyx(i)-F_Pyx(i);
        end
        
        % moddF/dT_P for pyroxenite
        if dFdT_P_Pyx(i)<0
            moddFdT_P_Pyx(i)=0.0001;
        else
            moddFdT_P_Pyx(i)=dFdT_P_Pyx(i);
        end
        
        % melt fraction of peridotite lithology
        F_Per(i)=FOT(i);
        
        % dFdT_P_Per peridotite
        dFdT_P_Per(i)=FTpp_Per(i)-F_Per(i);
        
        % dTdP_F peridotite
        dTdP_F_Per(i)=dTsdP_s_Per(i);
        if F_Per(i)<=0
            dTdP_F_PerT(i)=dTsdP_s_Per(i);
        else
            dTdP_F_PerT(i)=dTdP_F_Per(i);
        end
        
        % dTdP_F pyroxenite
        dTdP_F_Pyx(i)=dTsdP_s_Pyx(i);
        if F_Pyx(i)<=0
            dTdP_F_PyxT(i)=dTsdP_s_Pyx(i);
        else
            dTdP_F_PyxT(i)=dTdP_F_Pyx(i);
        end
        
        % pyroxenite extra values
        dFdP1_Pyx(i)=0;
        dFdP2_Pyx(i)=0;
        dFdP3_Pyx(i)=0;
        dFdP4_Pyx(i)=0;
        dFdP_S_Pyx(i)=0;
        
        % pyroxenite extra values
        dFdP1_Per(i)=0;
        dFdP2_Per(i)=0;
        dFdP3_Per(i)=0;
        dFdP_S_Per(i)=0;
        
        % change in temperature with pressure
        if (F_Per(i)==0) && (F_Pyx(i)==0)
            DTDP(i)=(T_SA(i+1)-T_SA(i))/(P(i+1)-P(i));
        else
            if F_Pyx(i)==0
                DTDP(i)=dTdP_F_PerT(i)+(1/dFdT_P_Per(i))*dFdP3_Per(i);
            else
                DTDP(i)=dTdP_F_PyxT(i)+(1/moddFdT_P_Pyx(i))*dFdP4_Pyx(i);
            end
        end
        
        % dT/dP_S
        if FTprimPyx(i)>=100
            dTdP_S_final(i)=dTdP_F_PerT(i)+(1/dFdT_P_Per(i))*dFdP3_Per(i);
        else
            dTdP_S_final(i)=DTDP(i);
        end
        
    else
        
        
        % temperature of solid mantle
        if (TsolPyx(i)>T_SA(i)) && (TsolPer(i)>T_SA(i))
            T(i)=T_SA(i);
        else
            T(i)=T(i-1)+dTdP_S_final(i-1)*(P(i)-P(i-1));
        end
        
        if T5Pyx(i)>Tcpx_outPyx(i) % calculating T' for the pyroxenite column
            TprimPyxT(i)=-1;
        else
            TprimPyxT(i)=(T(i)-T5Pyx(i))/(Tcpx_outPyx(i)-T5Pyx(i));
        end
        
        % calculating T' for the peridotite column
        TprimPerT(i)=(T(i)-TsolPer(i))/(Tliqlherz(i)-TsolPer(i));
        
        if TprimPyxT(i)<TprimPyx(i) % F(T') for pyroxenite
            FTprimPyx(i)=0;
        else
            if (((A+B*MgNo)*(P(i)^2)+(C+D*MgNo)*P(i)+E)*(TprimPyxT(i)^2)+(AA*CaO+BB)*TprimPyxT(i)+5)<0
                FTprimPyx(i)=0;
            else
                FTprimPyx(i)=(((A+B*MgNo)*(P(i)^2)+(C+D*MgNo)*P(i)+E)*(TprimPyxT(i)^2)+(AA*CaO+BB)*TprimPyxT(i)+5);
            end
        end
        
        if TprimPerT(i)<0 % F0-CPX (T') for peridotite
            F0_cpxT(i)=0;
        else
            F0_cpxT(i)=(TprimPerT(i)^(1.5))*100;
        end
        
        if T5Pyx(i)>Tcpx_outPyx(i) % T''(T+1) for pyroxenite
            Tpp1_Pyx(i)=-1;
        else
            Tpp1_Pyx(i)=TprimPyxT(i)+(1/(Tcpx_outPyx(i)-T5Pyx(i)));
        end
        
        % T''(T+1) for peridotite
        Tpp1_Per(i)=(T(i)+1-TsolPer(i))/(Tliqlherz(i)-TsolPer(i));
        
        % F-cpx(T'') peridotite
        if Tpp1_Per(i)<0
            F_cpxT(i)=0;
        else
            F_cpxT(i)=(Tpp1_Per(i)^(1.5))*100;
        end
        
        % F0(T') peridotite
        if F0_cpxT(i)<Fcpx_out(i)*100
            FOT(i)=F0_cpxT(i);
        else
            FOT(i)=(Fcpx_out(i)+((1-Fcpx_out(i))*((T(i)-Tcpx_outPer(i))/(TliqPer(i)-Tcpx_outPer(i))))^1.5)*100;
        end
        
        % F(T'') pyroxenite
        if Tpp1_Pyx(i)<TprimPyx(i)
            FTpp1_Pyx(i)=0;
        else
            if (((A+B*MgNo)*(P(i)^2)+(C+D*MgNo)*P(i)+E)*(Tpp1_Pyx(i)^2)+(AA*CaO+BB)*Tpp1_Pyx(i)+5)<0
                FTpp1_Pyx(i)=0;
            else
                FTpp1_Pyx(i)=(((A+B*MgNo)*(P(i)^2)+(C+D*MgNo)*P(i)+E)*(Tpp1_Pyx(i)^2)+(AA*CaO+BB)*Tpp1_Pyx(i)+5);
            end
        end
        
        % F(T'') peridotite
        if F_cpxT(i)<Fcpx_out(i)*100
            FTpp_Per(i)=F_cpxT(i);
        else
            FTpp_Per(i)=(Fcpx_out(i)+((1-Fcpx_out(i))*((T(i)+1-Tcpx_outPer(i))/(TliqPer(i)-Tcpx_outPer(i))))^1.5)*100;
        end
        
        if FTprimPyx(i-1)>0 % working out F_Pyx
            F_Pyx(i)=F_Pyx(i-1)+dFdP4_Pyx(i-1)*(P(i)-P(i-1));
        else
            F_Pyx(i)=FTprimPyx(i);
        end
        
        if FOT(i-1)<=0 % working out F_Per
            F_Per(i)=FOT(i);
        else
            F_Per(i)=F_Per(i-1)+dFdP3_Per(i-1)*(P(i)-P(i-1));
        end
        
        % dF/dT_P for pyroxenite (cpx exhaustion matters)
        if T(i)>Tcpx_outPyx(i)
            dFdT_P_Pyx(i)=0.3;
        else
            dFdT_P_Pyx(i)=FTpp1_Pyx(i)-F_Pyx(i);
        end
        
        % moddF/dT_P for pyroxenite
        if dFdT_P_Pyx(i)<0
            moddFdT_P_Pyx(i)=0.0001;
        else
            moddFdT_P_Pyx(i)=dFdT_P_Pyx(i);
        end
        
        % dFdT_P_Per peridotite
        dFdT_P_Per(i)=FTpp_Per(i)-F_Per(i);
        
                
        % dTdP_F peridotite - need to first work out Tm and T'm to
        % determine the dTdP_F_Per values
        TpM_Per(i-1)=(FOT(i)/100)^(1/1.5);
        if F_Per(i-1)<(Fcpx_out(i-1)*100)
            TM_Per(i-1)=TpM_Per(i-1)*(Tliqlherz(i-1)-TsolPer(i-1))+TsolPer(i-1);
        else
            TM_Per(i-1)=TpM_Per(i-1)*(TliqPer(i-1)-Tcpx_outPer(i-1))+Tcpx_outPer(i-1);
        end
        
        % now can calculate the dtdp_f parameter
        if TpM_Per(i-1)<TprimPerT(i-1)
            dTdP_F_Per(i)=dTdP_F_Per(i-1);
        else
            dTdP_F_Per(i)=(T(i)-TM_Per(i-1))/(P(i)-P(i-1));
        end
        
        % and therefore
        if F_Per(i)<=0
            dTdP_F_PerT(i)=dTsdP_s_Per(i);
        else
            dTdP_F_PerT(i)=dTdP_F_Per(i);
        end
        
        % now need to do the same for the pyroxenite component
        % calcl T'm and Tm parameters
        TpM_Pyx(i-1)=(-(AA*CaO+BB)+(((AA*CaO+BB)^2)-4*(5-FTprimPyx(i))*((A+B*MgNo)*(P(i-1)^2)+(C+D*MgNo)*P(i-1)+E))^0.5)/(2*((A+B*MgNo)*(P(i-1)^2)+(C+D*MgNo)*P(i-1)+E));
        TM_Pyx(i-1)=TpM_Pyx(i-1)*(Tcpx_outPyx(i-1)-T5Pyx(i-1))+T5Pyx(i-1);
        
        
        % dTdP_F pyroxenite - this is the issue!
        dTdP_F_Pyx(i)=dTsdP_s_Pyx(i-1);
        
        if F_Pyx(i)<=0
            dTdP_F_PyxT(i)=dTsdP_s_Pyx(i);
        else
            dTdP_F_PyxT(i)=dTdP_F_Pyx(i);
        end
        
        % now need to calculate the complicated shit for both peridotite
        % and pyroxenite and then finally calculate those two important
        % parameters again at the end...
        
        % starting with the peridotite shizzle
%         if (dFdT_P_Per(i)==0) || (moddFdT_P_Pyx(i)==0)
%             dFdP1_Per(i)=0;
%         else
            if (FTprimPyx(i)<=0) || (F_Pyx(i)>=100)
                dFdP1_Per(i)=((dTdP_F_PerT(i)-(alphaPd*(Tp+273.15)/(FracPyx*rhoPy*CpPy+(1-FracPyx)*rhoPd*CpPd)*10^9))/((1-FracPyx)*(Tp+273.15)*SPd/(FracPyx*CpPy+(1-FracPyx)*CpPd)+(100/dFdT_P_Per(i))))*(-100);
            else
                dFdP1_Per(i)=(dTdP_F_PerT(i)-((1-FracPyx)*alphaPd+FracPyx*alphaPy)*((Tp+273.15)/(FracPyx*rhoPy*CpPy+(1-FracPyx)*rhoPd*CpPd)*10^9)+FracPyx*(Tp+273.15)*(SPy/CpPy)*(dTdP_F_PerT(i)-dTdP_F_PyxT(i))/(100/moddFdT_P_Pyx(i)))/((1-FracPyx)*(Tp+273.15)*(SPd/CpPd)+FracPyx*(Tp+273.15)*(SPy/CpPy)*(100/dFdT_P_Per(i))/(100/moddFdT_P_Pyx(i))+(100/dFdT_P_Per(i)))*(-100);
            end
%         end
        
        if F_Per(i)<=0
            dFdP2_Per(i)=0;
        else
            dFdP2_Per(i)=dFdP1_Per(i);
        end
        
        if dFdT_P_Per(i)<0 || dFdP2_Per(i)>=0
            dFdP3_Per(i)=0;
        else
            dFdP3_Per(i)=dFdP2_Per(i);
        end
        
        dFdP_S_Per(i)=-dFdP3_Per(i);
        
        % and now for the pyroxenite
%         if (moddFdT_P_Pyx(i)==0) || (dFdT_P_Per(i)==0)
%             dFdP1_Pyx(i)=0;
%         else
            if dFdP_S_Per(i)<=0
                dFdP1_Pyx(i)=((dTdP_F_PyxT(i)-((alphaPy*(Tp+273.15)/(FracPyx*rhoPy*CpPy+(1-FracPyx)*rhoPd*CpPd))*10^9))/(FracPyx*(Tp+273.15)*SPy/(FracPyx*CpPy+(1-FracPyx)*CpPd)+(100/moddFdT_P_Pyx(i))))*(-100);
            else
                dFdP1_Pyx(i)=((dTdP_F_PyxT(i)-(FracPyx*alphaPy+(1-FracPyx)*alphaPd)*((Tp+273.15)/(FracPyx*rhoPy*CpPy+(1-FracPyx)*rhoPd*CpPd))*10^9+(1-FracPyx)*(Tp+273.15)*(SPd/CpPd)*(dTdP_F_PyxT(i)-dTdP_F_PerT(i))/(100/dFdT_P_Per(i)))/(FracPyx*(Tp+273.15)*(SPy/CpPy)+(1-FracPyx)*(Tp+273.15)*(SPd/CpPd)*(100/moddFdT_P_Pyx(i))/(100/dFdT_P_Per(i))+(100/moddFdT_P_Pyx(i))))*(-100);
            end
%         end
        
        if moddFdT_P_Pyx(i)==0
            dFdP2_Pyx(i)=0;
        else 
            if FPyx_SA(i+1)==0
                dFdP2_Pyx(i)=0;
            else
                dFdP2_Pyx(i)=dFdP1_Pyx(i);
            end
        end
        
        if dFdP2_Pyx(i)>=0
            dFdP3_Pyx(i)=0;
        else
            dFdP3_Pyx(i)=dFdP2_Pyx(i);
        end
        
        if (F_Pyx(i)>100)
            dFdP4_Pyx(i)=0;
        else
            dFdP4_Pyx(i)=dFdP3_Pyx(i);
        end
        
        dFdP_S_Pyx(i)=-dFdP4_Pyx(i);
        
        % calculate the final two variables
        % change in temperature with pressure
        if (F_Per(i)==0) && (F_Pyx(i)==0)
            DTDP(i)=(T_SA(i+1)-T_SA(i))/(P(i+1)-P(i));
        else
            if F_Pyx(i)==0
                DTDP(i)=dTdP_F_PerT(i)+(1/dFdT_P_Per(i))*dFdP3_Per(i);
            else
                DTDP(i)=dTdP_F_PyxT(i)+(1/moddFdT_P_Pyx(i))*dFdP4_Pyx(i);
            end
        end
        
        % dT/dP_S
        if FTprimPyx(i)>=100
            dTdP_S_final(i)=dTdP_F_PerT(i)+(1/dFdT_P_Per(i))*dFdP3_Per(i);
        else
            dTdP_S_final(i)=DTDP(i);
        end
    end
end
        
Pressure=P(1:end-1);
T_solidadiabat=T_SA(1:end-1);

end