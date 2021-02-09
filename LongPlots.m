%% Code to run meltPx function and add in code to calculate melt chemistry
clear all
close all

%% input parameters
% Pyroxenite composition
SiO2=47.5395553422756; TiO2=0.3832220670177; Al2O3=7.14067512672454;
Cr2O3=0.0464700082784867; FeO=7.38957897155817; MnO=0.165943726786654;
MgO=24.3001432174423; CaO=11.4300141859886; Na2O=1.56690635075243;
K2O=0.0374910031755343;

dFePdstart=0.01; dFePystart=0.20; % required for the functions to run fully.

N=500; % extent of mixing for Dirichlet function
m=0; % select 1 for channelised flow, 0 for no channelised flow

% load data
Data=readtable('Data/GSCCorrectedFinal.xlsx');

% user input parameters for trace code
Dsource=3; % 1=WHDD 2% melt extraction - 2=WHDD+Donelly - 3=WHDD depleted
Esource=2; % 1-SMPM - 2=PW09+WHDD - 3=KG1 Lambart 2017

% normalisation
nRb = 0.635; nBa = 6.989; nTh = 0.085; nU = 0.021; nNb = 0.713; nTa = 0.041; nLa = 0.687; nCerium = 1.775; nPb = 0.071; nPr = 0.276; nSr = 21.1; nNd = 1.354; 
nZr = 11.2; nHf = 0.309; nSm = 0.444; nEu = 0.168; nGd = 0.596; nTb = 0.108; nDy = 0.737; nHo = 0.164; nEr = 0.48; nYb = 0.493; nY = 4.55; nLu = 0.074;
nH2O = 350; nFlo=30;

norm=[nLa nCerium nPr nNd nSm nEu nGd nTb nDy nHo nEr nYb nLu];

if m==0

%% If doing longitude variations
% EGSC
for uuu=1:200
Longplot(uuu)=(uuu-1)*(90.8-86)/200+86;
Ur=exp(-((90.8-Longplot(uuu))*0.75))*10; %10
Ursave(uuu)=Ur+1;
Tp=(12/200)*(uuu-1)+1347; %exp(-((90.8-Longplot(uuu))*0.05))*40+1310;
Tpsave(uuu)=Tp;

% user input parameters
FracPyx=0.03; % Fraction of pyroxenite in the mantle source
FracCpx=0.15; % Mass fraction of cpx in the subsolidus peridotite
    
dtop=0.3; % base of lithosphere/top of melt column
Dchange=2; % prssure (GPa) at which melt region changes from triangular to rectangular (melt extraction region)
Ffactor=0.75; % extent of melt remaining following fractional crystallisation of pure mantle melts
m=0; % mass fraction of channelised melt
Xd=1-FracPyx;

% MELT_PX
[Pressure,F_Per,F_Pyx,T,T_solidadiabat,TsolPer,TsolPyx]=MeltPX(Tp,FracPyx,FracCpx,SiO2,TiO2,Al2O3,Cr2O3,FeO,MnO,MgO,CaO,Na2O,K2O);
% Uchange=Pressure(find(F_Per>0,1)); % depth at which changes in the relative upwelling velocity starts to occur
% Uchange=1.2; %(round(Uchange.*100))./100;
Uchange=Pressure(find(F_Per>0,1)); % depth at which changes in the relative upwelling velocity starts to occur
Uchange=(round(Uchange.*100))./100;

% melt PX calculation code
Dsource=3;
[Ctfinal,Pdd,CREEf,CREE,CfH2O,Hflux,Urin,PPP,Cmix]=MeltPXtraceDirichletFe(Pressure,F_Pyx,F_Per,T,Dsource,Esource,dtop,Ur,Uchange,Dchange,Xd,Ffactor,m,N,dFePdstart,dFePystart);

LaYb(uuu)=(CREE(1)/nLa)/(CREE(end-1)/nYb);
CTsave(uuu)=Ctfinal;
meanLaYb(uuu)=mean((Cmix(:,3)./nLa)./(Cmix(:,end-3)./nYb));
stdLaYb(uuu)=std((Cmix(:,3)./nLa)./(Cmix(:,end-3)./nYb));
HF(uuu)=Hflux;
meanHF(uuu)=mean((Cmix(:,end-1)*10^(-6))*0.057*Ctfinal*2900);
stdHF(uuu)=std((Cmix(:,end-1)*10^(-6))*0.057*Ctfinal*2900);
meanH2O(uuu)=mean(Cmix(:,end-1)./Ffactor);
stdH2O(uuu)=std(Cmix(:,end-1)./Ffactor);
end  

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse','YAxisLocation','right')
plot(Longplot,LaYb,'-k')
plot(Longplot,meanLaYb,'-k','LineWidth',2)
plot(Longplot,meanLaYb+2*stdLaYb,'-.r','LineWidth',2)
plot(Longplot,meanLaYb-2*stdLaYb,'-.r','LineWidth',2)
xlabel('Longitude (^{o}W)','FontSize',16)
ylabel('La_{n}/Yb_{n}','FontSize',16)
box on
xlim([85 91])
for i=1:length(Data.La)
    if (Data.Long(i)>(-90.8)) && (Data.Long(i)<-85.5)
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)<0.8
        plot((-1)*(Data.Long(i)),(Data.La(i)/nLa)/(Data.Yb(i)/nYb),'ok','MarkerFaceColor',[1 1 1],'MarkerSize',8)
    end
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)>1.2
        plot((-1)*(Data.Long(i)),(Data.La(i)/nLa)/(Data.Yb(i)/nYb),'ok','MarkerFaceColor','green','MarkerSize',8)
    end
    if ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))>0.8 && ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))<1.2
        plot((-1)*(Data.Long(i)),(Data.La(i)/nLa)/(Data.Yb(i)/nYb),'ok','MarkerFaceColor','blue','MarkerSize',8)
    end
    end
end

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse','YAxisLocation','right')
plot(Longplot,(CfH2O./Ffactor)./10000,'-k')
plot(Longplot,(meanH2O)./10000,'-k','LineWidth',2)
plot(Longplot,(meanH2O+2*stdH2O)./10000,'-.r','LineWidth',2)
plot(Longplot,(meanH2O-2*stdH2O)./10000,'-.r','LineWidth',2)
xlabel('Longitude (^{o}W)','FontSize',16)
ylabel('H_{2}O_{[8]}^{*} (wt%)','FontSize',16)
box on
ylim([0 0.8])
xlim([85 91])
for i=1:length(Data.La)
    if (Data.Long(i)>(-90.8)) && (Data.Long(i)<-85.5)
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)<0.8
        plot((-1)*(Data.Long(i)),(Data.H2Ofull(i)),'ok','MarkerFaceColor',[1 1 1],'MarkerSize',8)
    end
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)>1.2
        plot((-1)*(Data.Long(i)),(Data.H2Ofull(i)),'ok','MarkerFaceColor','green','MarkerSize',8)
    end
    if ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))>0.8 && ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))<1.2
        plot((-1)*(Data.Long(i)),(Data.H2Ofull(i)),'ok','MarkerFaceColor','blue','MarkerSize',8)
    end
    end
end

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse','YAxisLocation','right')
plot(Longplot,HF,'-k')
plot(Longplot,meanHF,'-k','LineWidth',2)
plot(Longplot,meanHF+2*stdHF,'-.r','LineWidth',2)
plot(Longplot,meanHF-2*stdHF,'-.r','LineWidth',2)
xlabel('Longitude (^{o}W)','FontSize',16)
ylabel('H_{2}O flux (kg.m^{-1}.yr^{-1})','FontSize',16)
box on
xlim([85 91])

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse','YAxisLocation','right')
plot(Longplot,CTsave,'-k')
xlabel('Longitude (^{o}W)','FontSize',16)
box on
xlim([85 91])
ylabel('Crustal Thickness','FontSize',16)
ylim([6500 11000])

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse','YAxisLocation','right')
plot(Longplot,Tpsave,'-k')
xlabel('Longitude (^{o}W)','FontSize',16)
box on
xlim([85 91])
ylabel('Mantle Potential Temperature (^{o}C)','FontSize',16)

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse','YAxisLocation','right')
plot(Longplot,Ursave,'-k')
xlabel('Longitude (^{o}W)','FontSize',16)
box on
xlim([85 91])
ylabel('U_{r-max}','FontSize',16)
ylim([0 12])

% WGSC
for uuu=1:200
Longplot(uuu)=(uuu-1)*(90.8-95.5)/200+95.5;
Ur=exp(-((Longplot(uuu)-90.8)*0.6))*5;
Ursave(uuu)=Ur+1;
Tp=(20/200)*(uuu-1)+1337; %exp(-((90.8-Longplot(uuu))*0.05))*40+1310;
Tpsave(uuu)=Tp;

% user input parameters
FracPyx=0.03; % Fraction of pyroxenite in the mantle source
FracCpx=0.15; % Mass fraction of cpx in the subsolidus peridotite
    
dtop=0.3; % base of lithosphere/top of melt column
Dchange=2; % prssure (GPa) at which melt region changes from triangular to rectangular (melt extraction region)
Ffactor=0.75; % extent of melt remaining following fractional crystallisation of pure mantle melts
m=0; % mass fraction of channelised melt
Xd=1-FracPyx;

% MELT_PX
[Pressure,F_Per,F_Pyx,T,T_solidadiabat,TsolPer,TsolPyx]=MeltPX(Tp,FracPyx,FracCpx,SiO2,TiO2,Al2O3,Cr2O3,FeO,MnO,MgO,CaO,Na2O,K2O);
Uchange=Pressure(find(F_Per>0,1)); % depth at which changes in the relative upwelling velocity starts to occur
Uchange=(round(Uchange.*100))./100;

% melt PX calculation code
Dsource=4;
[Ctfinal,Pdd,CREEf,CREE,CfH2O,Hflux,Urin,PPP,Cmix]=MeltPXtraceDirichletFe(Pressure,F_Pyx,F_Per,T,Dsource,Esource,dtop,Ur,Uchange,Dchange,Xd,Ffactor,m,N,dFePdstart,dFePystart);

LaYb(uuu)=(CREE(1)/nLa)/(CREE(end-1)/nYb);
CTsave(uuu)=Ctfinal;
meanLaYb(uuu)=mean((Cmix(:,3)./nLa)./(Cmix(:,end-3)./nYb));
stdLaYb(uuu)=std((Cmix(:,3)./nLa)./(Cmix(:,end-3)./nYb));
HF(uuu)=Hflux;
meanHF(uuu)=mean((Cmix(:,end-1)*10^(-6))*0.057*Ctfinal*2900);
stdHF(uuu)=std((Cmix(:,end-1)*10^(-6))*0.057*Ctfinal*2900);
meanH2O(uuu)=mean(Cmix(:,end-1)./Ffactor);
stdH2O(uuu)=std(Cmix(:,end-1)./Ffactor);
end  

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse')
plot(Longplot,LaYb,'-k')
plot(Longplot,meanLaYb,'-k','LineWidth',2)
plot(Longplot,meanLaYb+2*stdLaYb,'-.r','LineWidth',2)
plot(Longplot,meanLaYb-2*stdLaYb,'-.r','LineWidth',2)
xlabel('Longitude (^{o}W)','FontSize',16)
ylabel('La_{n}/Yb_{n}','FontSize',16)
box on
xlim([90 96])
for i=1:length(Data.La)
    if (Data.Long(i)<(-90.8)) && (Data.Long(i)>-95.5)
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)<0.8
        plot((-1)*(Data.Long(i)),(Data.La(i)/nLa)/(Data.Yb(i)/nYb),'ok','MarkerFaceColor',[1 1 1],'MarkerSize',8)
    end
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)>1.2
        plot((-1)*(Data.Long(i)),(Data.La(i)/nLa)/(Data.Yb(i)/nYb),'ok','MarkerFaceColor','green','MarkerSize',8)
    end
    if ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))>0.8 && ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))<1.2
        plot((-1)*(Data.Long(i)),(Data.La(i)/nLa)/(Data.Yb(i)/nYb),'ok','MarkerFaceColor','blue','MarkerSize',8)
    end
    end
end

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse')
plot(Longplot,(CfH2O./Ffactor)./10000,'-k')
plot(Longplot,(meanH2O)./10000,'-k','LineWidth',2)
plot(Longplot,(meanH2O+2*stdH2O)./10000,'-.r','LineWidth',2)
plot(Longplot,(meanH2O-2*stdH2O)./10000,'-.r','LineWidth',2)
xlabel('Longitude (^{o}W)','FontSize',16)
ylabel('H_{2}O_{[8]}^{*} (wt%)','FontSize',16)
box on
xlim([90 96])
ylim([0 0.8])
for i=1:length(Data.La)
    if (Data.Long(i)<(-90.8)) && (Data.Long(i)>-95.5)
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)<0.8
        plot((-1)*(Data.Long(i)),(Data.H2Ofull(i)),'ok','MarkerFaceColor',[1 1 1],'MarkerSize',8)
    end
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)>1.2
        plot((-1)*(Data.Long(i)),(Data.H2Ofull(i)),'ok','MarkerFaceColor','green','MarkerSize',8)
    end
    if ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))>0.8 && ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))<1.2
        plot((-1)*(Data.Long(i)),(Data.H2Ofull(i)),'ok','MarkerFaceColor','blue','MarkerSize',8)
    end
    end
end

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse')
plot(Longplot,HF,'-k')
plot(Longplot,meanHF,'-k','LineWidth',2)
plot(Longplot,meanHF+2*stdHF,'-.r','LineWidth',2)
plot(Longplot,meanHF-2*stdHF,'-.r','LineWidth',2)
xlabel('Longitude (^{o}W)','FontSize',16)
ylabel('H_{2}O flux (kg.m^{-1}.yr^{-1})','FontSize',16)
box on
xlim([90 96])

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse')
plot(Longplot,CTsave,'-k')
xlabel('Longitude (^{o}W)','FontSize',16)
box on
xlim([90 96])
ylabel('Crustal Thickness','FontSize',16)
ylim([6500 11000])

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse')
plot(Longplot,Tpsave,'-k')
xlabel('Longitude (^{o}W)','FontSize',16)
box on
xlim([90 96])
ylabel('Mantle Potential Temperature (^{o}C)','FontSize',16)

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse')
plot(Longplot,Ursave,'-k')
xlabel('Longitude (^{o}W)','FontSize',16)
box on
xlim([90 96])
ylabel('U_{r-max}','FontSize',16)
ylim([0 12])
end

if m>0
for uuu=1:200
Longplot(uuu)=(uuu-1)*(90.8-86)/200+86;
% Ur=0;
Ur=exp(-((90.8-Longplot(uuu))*0.6))*2.2;
Ursave(uuu)=Ur+1;
m=exp(-((90.8-Longplot(uuu))*0.9))*0.16-0.002;
msave(uuu)=m;
Tp=-0.6906*(Longplot(uuu)^2)+124.67*Longplot(uuu)-4269.8; %-1.3908*(Longplot(uuu)^2)+250.84*Longplot(uuu)-9945.2; %-0.497*(Longplot(uuu)^2)+89.925*Longplot(uuu)-2715.7;%-0.6906*(Longplot(uuu)^2)+124.67*Longplot(uuu)-4271.6; %(12/200)*(uuu-1)+1342; %exp(-((90.8-Longplot(uuu))*0.05))*40+1310;
Tpsave(uuu)=Tp;

% user input parameters
FracPyx=0.03; % Fraction of pyroxenite in the mantle source
FracCpx=0.15; % Mass fraction of cpx in the subsolidus peridotite
    
dtop=exp(-((90.8-Longplot(uuu))*2.2))*0.12+0.3;
dtop=round(dtop*100)/100;% base of lithosphere/top of melt column
Dchange=2.2; % prssure (GPa) at which melt region changes from triangular to rectangular (melt extraction region)
Ffactor=0.74; % extent of melt remaining following fractional crystallisation of pure mantle melts
Xd=1-FracPyx;

% MELT_PX
[Pressure,F_Per,F_Pyx,T,T_solidadiabat,TsolPer,TsolPyx]=MeltPX(Tp,FracPyx,FracCpx,SiO2,TiO2,Al2O3,Cr2O3,FeO,MnO,MgO,CaO,Na2O,K2O);
Uchange=Pressure(find(F_Per>0,1)); % depth at which changes in the relative upwelling velocity starts to occur
Uchange=(round(Uchange.*100))./100;

% melt PX calculation code
Dsource=3;
[Ctfinal,Pdd,CREEf,CREE,CfH2O,Hflux,Urin,PPP,Cmix,CH2OGarnet]=MeltPXtraceDirichletFe(Pressure,F_Pyx,F_Per,T,Dsource,Esource,dtop,Ur,Uchange,Dchange,Xd,Ffactor,m,N,dFePdstart,dFePystart);

LaYb(uuu)=(CREE(1)/nLa)/(CREE(end-1)/nYb);
CTsave(uuu)=Ctfinal;
meanLaYb(uuu)=mean((Cmix(:,3)./nLa)./(Cmix(:,end-3)./nYb));
stdLaYb(uuu)=std((Cmix(:,3)./nLa)./(Cmix(:,end-3)./nYb));
HF(uuu)=Hflux;
meanHF(uuu)=mean((Cmix(:,end-1)*10^(-6))*0.057*Ctfinal*2900);
stdHF(uuu)=std((Cmix(:,end-1)*10^(-6))*0.057*Ctfinal*2900);
meanH2O(uuu)=mean(Cmix(:,end-1)./Ffactor);
stdH2O(uuu)=std(Cmix(:,end-1)./Ffactor);
PropH2O(uuu)=(m*CH2OGarnet)/CfH2O;
end  

meanPropH2OEGSC=mean(PropH2O);

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse','YAxisLocation','right')
plot(Longplot,LaYb,'-k')
plot(Longplot,meanLaYb,'-k','LineWidth',2)
plot(Longplot,meanLaYb+2*stdLaYb,'-.r','LineWidth',2)
plot(Longplot,meanLaYb-2*stdLaYb,'-.r','LineWidth',2)
xlabel('Longitude (^{o}W)','FontSize',16)
ylabel('La_{n}/Yb_{n}','FontSize',16)
box on
xlim([85 91])
for i=1:length(Data.La)
    if (Data.Long(i)>(-90.8)) && (Data.Long(i)<-85.5)
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)<0.8
        plot((-1)*(Data.Long(i)),(Data.La(i)/nLa)/(Data.Yb(i)/nYb),'ok','MarkerFaceColor',[1 1 1],'MarkerSize',8)
    end
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)>1.2
        plot((-1)*(Data.Long(i)),(Data.La(i)/nLa)/(Data.Yb(i)/nYb),'ok','MarkerFaceColor','green','MarkerSize',8)
    end
    if ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))>0.8 && ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))<1.2
        plot((-1)*(Data.Long(i)),(Data.La(i)/nLa)/(Data.Yb(i)/nYb),'ok','MarkerFaceColor','blue','MarkerSize',8)
    end
    end
end

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse','YAxisLocation','right')
plot(Longplot,(CfH2O./Ffactor)./10000,'-k')
plot(Longplot,(meanH2O)./10000,'-k','LineWidth',2)
plot(Longplot,(meanH2O+2*stdH2O)./10000,'-.r','LineWidth',2)
plot(Longplot,(meanH2O-2*stdH2O)./10000,'-.r','LineWidth',2)
xlabel('Longitude (^{o}W)','FontSize',16)
ylabel('H_{2}O_{[8]}^{*} (wt%)','FontSize',16)
box on
ylim([0 0.8])
xlim([85 91])
for i=1:length(Data.La)
    if (Data.Long(i)>(-90.8)) && (Data.Long(i)<-85.5)
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)<0.8
        plot((-1)*(Data.Long(i)),(Data.H2Ofull(i)),'ok','MarkerFaceColor',[1 1 1],'MarkerSize',8)
    end
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)>1.2
        plot((-1)*(Data.Long(i)),(Data.H2Ofull(i)),'ok','MarkerFaceColor','green','MarkerSize',8)
    end
    if ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))>0.8 && ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))<1.2
        plot((-1)*(Data.Long(i)),(Data.H2Ofull(i)),'ok','MarkerFaceColor','blue','MarkerSize',8)
    end
    end
end

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse','YAxisLocation','right')
plot(Longplot,HF,'-k')
plot(Longplot,meanHF,'-k','LineWidth',2)
plot(Longplot,meanHF+2*stdHF,'-.r','LineWidth',2)
plot(Longplot,meanHF-2*stdHF,'-.r','LineWidth',2)
xlabel('Longitude (^{o}W)','FontSize',16)
ylabel('H_{2}O flux (kg.m^{-1}.yr^{-1})','FontSize',16)
box on
xlim([85 91])

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse','YAxisLocation','right')
plot(Longplot,CTsave,'-k')
xlabel('Longitude (^{o}W)','FontSize',16)
box on
xlim([85 91])
ylabel('Crustal Thickness','FontSize',16)
ylim([6500 11000])

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse','YAxisLocation','right')
plot(Longplot,PropH2O,'-k')
xlabel('Longitude (^{o}W)','FontSize',16)
box on
xlim([85 91])
ylabel('Proportion of H_{2}O from channelised melts','FontSize',16)

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse','YAxisLocation','right')
plot(Longplot,msave,'-k')
xlabel('Longitude (^{o}W)','FontSize',16)
box on
xlim([85 91])
ylabel('Mass of channelised melt','FontSize',16)
ylim([0 0.2])

% % WGSC
Longplot=linspace(90.8,95.5,200);
Y=(1./(0.12.*sqrt(2.*(pi))).*exp((-(Longplot-92.25).^2)./(2.*(0.12.^2))));
Y2=(1./(0.12.*sqrt(2.*(pi))).*exp((-(Longplot-91.8).^2)./(2.*(0.12.^2))));
Y3=(1./(0.12.*sqrt(2.*(pi))).*exp((-(Longplot-91.3).^2)./(2.*(0.12.^2))));
Ysum=Y+Y2+0.8.*Y3;
msave=(Ysum./max(Ysum)).*0.12;
Tpsave=-0.4318.*(Longplot.^2)+77.676.*Longplot-2147.0; %1330 to 1345

for uuu=1:200
Ur=exp(-((Longplot(uuu)-90.8)*0.5))*1.3;
Ursave(uuu)=Ur+1;
% Longplot(uuu)=(uuu-1)*(91.8-95.5)/200+95.5;
% Ur=exp(-((Longplot(uuu)-91.8)*2))*6;
% Ursave(uuu)=Ur+1;
m=msave(uuu);

Tp=Tpsave(uuu);
% Tp=-(15/200)*(uuu-1)+1347; %exp(-((90.8-Longplot(uuu))*0.05))*40+1310;
% Tpsave(uuu)=Tp;

% user input parameters
FracPyx=0.035; % Fraction of pyroxenite in the mantle source0.
FracCpx=0.15; % Mass fraction of cpx in the subsolidus peridotite
    
dtop=0.3; % base of lithosphere/top of melt column
Dchange=2; % prssure (GPa) at which melt region changes from triangular to rectangular (melt extraction region)
Ffactor=0.75; % extent of melt remaining following fractional crystallisation of pure mantle melts

Xd=1-FracPyx;

% MELT_PX
[Pressure,F_Per,F_Pyx,T,T_solidadiabat,TsolPer,TsolPyx]=MeltPX(Tp,FracPyx,FracCpx,SiO2,TiO2,Al2O3,Cr2O3,FeO,MnO,MgO,CaO,Na2O,K2O);
Uchange=Pressure(find(F_Per>0,1)); % depth at which changes in the relative upwelling velocity starts to occur
Uchange=(round(Uchange.*100))./100;

% melt PX calculation code
Dsource=4;
[Ctfinal,Pdd,CREEf,CREE,CfH2O,Hflux,Urin,PPP,Cmix,CH2OGarnet]=MeltPXtraceDirichletFe(Pressure,F_Pyx,F_Per,T,Dsource,Esource,dtop,Ur,Uchange,Dchange,Xd,Ffactor,m,N,dFePdstart,dFePystart);

LaYb(uuu)=(CREE(1)/nLa)/(CREE(end-1)/nYb);
CTsave(uuu)=Ctfinal;
meanLaYb(uuu)=mean((Cmix(:,3)./nLa)./(Cmix(:,end-3)./nYb));
stdLaYb(uuu)=std((Cmix(:,3)./nLa)./(Cmix(:,end-3)./nYb));
HF(uuu)=Hflux;
meanHF(uuu)=mean((Cmix(:,end-1)*10^(-6))*0.057*Ctfinal*2900);
stdHF(uuu)=std((Cmix(:,end-1)*10^(-6))*0.057*Ctfinal*2900);
meanH2O(uuu)=mean(Cmix(:,end-1)./Ffactor);
stdH2O(uuu)=std(Cmix(:,end-1)./Ffactor);
PropH2O(uuu)=(m*CH2OGarnet)/CfH2O;
end  

% calculate proportion of H2O from channelised melts between 92.5 and
% 90.8oW
meanPropH2O=mean(PropH2O(1:73));

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse')
plot(Longplot,LaYb,'-k')
plot(Longplot,meanLaYb,'-k','LineWidth',2)
plot(Longplot,meanLaYb+2*stdLaYb,'-.r','LineWidth',2)
plot(Longplot,meanLaYb-2*stdLaYb,'-.r','LineWidth',2)
xlabel('Longitude (^{o}W)','FontSize',16)
ylabel('La_{n}/Yb_{n}','FontSize',16)
box on
xlim([90 96])
for i=1:length(Data.La)
    if (Data.Long(i)<(-90.8)) && (Data.Long(i)>-95.5)
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)<0.8
        plot((-1)*(Data.Long(i)),(Data.La(i)/nLa)/(Data.Yb(i)/nYb),'ok','MarkerFaceColor',[1 1 1],'MarkerSize',8)
    end
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)>1.2
        plot((-1)*(Data.Long(i)),(Data.La(i)/nLa)/(Data.Yb(i)/nYb),'ok','MarkerFaceColor','green','MarkerSize',8)
    end
    if ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))>0.8 && ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))<1.2
        plot((-1)*(Data.Long(i)),(Data.La(i)/nLa)/(Data.Yb(i)/nYb),'ok','MarkerFaceColor','blue','MarkerSize',8)
    end
    end
end

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse')
plot(Longplot,(CfH2O./Ffactor)./10000,'-k')
plot(Longplot,(meanH2O)./10000,'-k','LineWidth',2)
plot(Longplot,(meanH2O+2*stdH2O)./10000,'-.r','LineWidth',2)
plot(Longplot,(meanH2O-2*stdH2O)./10000,'-.r','LineWidth',2)
xlabel('Longitude (^{o}W)','FontSize',16)
ylabel('H_{2}O_{[8]}^{*} (wt%)','FontSize',16)
box on
xlim([90 96])
ylim([0 0.8])
for i=1:length(Data.La)
    if (Data.Long(i)<(-90.8)) && (Data.Long(i)>-95.5)
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)<0.8
        plot((-1)*(Data.Long(i)),(Data.H2Ofull(i)),'ok','MarkerFaceColor',[1 1 1],'MarkerSize',8)
    end
    if (Data.La(i)/nLa)/(Data.Sm(i)/nSm)>1.2
        plot((-1)*(Data.Long(i)),(Data.H2Ofull(i)),'ok','MarkerFaceColor','green','MarkerSize',8)
    end
    if ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))>0.8 && ((Data.La(i)/nLa)/(Data.Sm(i)/nSm))<1.2
        plot((-1)*(Data.Long(i)),(Data.H2Ofull(i)),'ok','MarkerFaceColor','blue','MarkerSize',8)
    end
    end
end

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse')
plot(Longplot,HF,'-k')
plot(Longplot,meanHF,'-k','LineWidth',2)
plot(Longplot,meanHF+2*stdHF,'-.r','LineWidth',2)
plot(Longplot,meanHF-2*stdHF,'-.r','LineWidth',2)
xlabel('Longitude (^{o}W)','FontSize',16)
ylabel('H_{2}O flux (kg.m^{-1}.yr^{-1})','FontSize',16)
box on
xlim([90 96])
ylim([0 5000])

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse')
plot(Longplot,CTsave,'-k')
xlabel('Longitude (^{o}W)','FontSize',16)
box on
xlim([90 96])
ylabel('Crustal Thickness','FontSize',16)
ylim([6500 11000])

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse')
plot(Longplot,PropH2O,'-k')
xlabel('Longitude (^{o}W)','FontSize',16)
box on
xlim([90 96])
ylabel('X_H{2}O from channelised melts','FontSize',16)

figure('rend','painters','pos',[10 10 550 300])
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'XDir','reverse')
plot(Longplot,msave,'-k')
xlabel('Longitude (^{o}W)','FontSize',16)
box on
xlim([90 96])
ylabel('Mass of channelised melt','FontSize',16)
ylim([0 0.2])    
    
end



