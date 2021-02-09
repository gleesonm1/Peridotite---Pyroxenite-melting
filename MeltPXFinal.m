%% Code to run meltPx function and add in code to calculate melt chemistry. 
% To run this code also hae dirichletRnd.m, MeltPX.m and
% MeltPXtraceDirichletFe.m in the same folder. There are some redundant
% features of these codes that are not used in detail in this study. The
% composition of samples 6D-1g and 6D-2g are included in the compressed
% (zipped) folder and loaded into this code to give an example of how this
% is used.

clear all
close all

%% input parameters for all codes, please follow instructions carefully
% Pyroxenite composition
SiO2=47.5395553422756; TiO2=0.3832220670177; Al2O3=7.14067512672454;
Cr2O3=0.0464700082784867; FeO=7.38957897155817; MnO=0.165943726786654;
MgO=24.3001432174423; CaO=11.4300141859886; Na2O=1.56690635075243;
K2O=0.0374910031755343;
% SiO2=44.1020; TiO2=0.5710; Al2O3=8.1878;
% Cr2O3=0.0280; FeO=12.2802; MnO=0.2571;
% MgO=18.888; CaO=14.9307; Na2O=0.7358;
% K2O=0.0194;

dFePdstart=0.01; dFePystart=0.20;

N=500; % extent of mixing for Dirichlet function
m=0; % select 1 for channelised flow, 0 for no channelised flow

% load data
Data=readtable('Data/6Donly.xlsx');

% user input parameters for trace code
NumberOfIterations=5000;
Dsource=4; % 1=WHDD 2% melt extraction - 2=WHDD+Donelly - 3=WHDD depleted
Esource=2; % 1-SMPM - 2=PW09+WHDD - 3=KG1 Lambart 2017

% normalisation
nRb = 0.635; nBa = 6.989; nTh = 0.085; nU = 0.021; nNb = 0.713; nTa = 0.041; nLa = 0.687; nCerium = 1.775; nPb = 0.071; nPr = 0.276; nSr = 21.1; nNd = 1.354; 
nZr = 11.2; nHf = 0.309; nSm = 0.444; nEu = 0.168; nGd = 0.596; nTb = 0.108; nDy = 0.737; nHo = 0.164; nEr = 0.48; nYb = 0.493; nY = 4.55; nLu = 0.074;
nH2O = 350; nFlo=30;

norm=[nLa nCerium nPr nNd nSm nEu nGd nTb nDy nHo nEr nYb nLu];

% REE of selected data
ResEE=[Data.La Data.Ce Data.Pr Data.Nd Data.Sm Data.Eu Data.Gd Data.Tb Data.Dy Data.Ho Data.Er Data.Yb Data.Lu];

% average and mean
for i=1:size(ResEE,2)
ResEa(i)=mean(ResEE(:,i));
ResEb(i)=0.2*mean(ResEE(:,i));%2.*std(ResEE(:,i));
end

figure('rend','painters','pos',[10 10 1200 400])
subaxis(1,3,1,'SpacingVert',0.04,'SpacingHoriz',0.04)
semilogy(ResEa./norm,'->k','MarkerFaceColor','green','MarkerSize',5,'MarkerSize',8)
hold on
xticklabels({'La','Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Yb','Lu'})
xticks([1 2 3 4 5 6 7 8 9 10 11 12 13])
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'YAxisLocation','left')
ylabel('Primitive mantle normalised','FontSize',16)
ylim([2 30])

% running the iterations
for uuu=1:NumberOfIterations
% user input parameters
Tp=(1365-1335)*rand(1,1)+1335; % oC mantle potential temperature
FracPyx=(0.05-0.03)*rand(1,1)+0.03; % Fraction of pyroxenite in the mantle source
FracCpx=0.15; % Mass fraction of cpx in the subsolidus peridotite
    
dtop=(0.8-0.3)*rand(1,1)+0.3; % base of lithosphere/top of melt column
dtop=(round(dtop.*100))./100;
Ur=40*rand(1,1); % relative rate of upwelling (1+Ur)
Dchange=2; % prssure (GPa) at which melt region changes from triangular to rectangular (melt extraction region)
Ffactor=(0.8-0.70)*rand(1,1)+0.70; % extent of melt remaining following fractional crystallisation of pure mantle melts
m=0; % mass fraction of channelised melt
Xd=1-FracPyx;

% MELT_PX
[Pressure,F_Per,F_Pyx,T,T_solidadiabat,TsolPer,TsolPyx]=MeltPX(Tp,FracPyx,FracCpx,SiO2,TiO2,Al2O3,Cr2O3,FeO,MnO,MgO,CaO,Na2O,K2O);

Uchange=(Pressure(find(F_Per>0,1))-(Pressure(find(F_Per>0,1))-0.1))*rand(1,1)+(Pressure(find(F_Per>0,1))-0.1); % depth at which changes in the relative upwelling velocity starts to occur
Uchange=(round(Uchange.*100))./100;

% melt PX calculation code
[Ctfinal,Pdd,CREEf,CREE,CfH2O,Hflux,Urin,PPP,Cmix,CH2OGarnet,Feiso]=MeltPXtraceDirichletFe(Pressure,F_Pyx,F_Per,T,Dsource,Esource,dtop,Ur,Uchange,Dchange,Xd,Ffactor,m,N,dFePdstart,dFePystart);


GF=(sum((abs(CREEf-ResEa))./(ResEb)));
Goodness=exp((-GF));

Penriched(uuu,1)=(Pdd);
dtopFINAL(uuu,1)=dtop;
FfactorFINAL(uuu,1)=Ffactor;
UrinFINAL(uuu,1)=Ur+1;
UchangeFINAL(uuu,1)=Uchange;
XdFINAL(uuu,1)=Xd;
CTFINAL(uuu,1)=Ctfinal;
HFINAL(uuu,1)=Hflux;
TpFINAL(uuu,1)=Tp;
FeFINAL(uuu,1)=Feiso;

Penriched(uuu,2)=Goodness;
dtopFINAL(uuu,2)=Goodness;
FfactorFINAL(uuu,2)=Goodness;
UrinFINAL(uuu,2)=Goodness;
UchangeFINAL(uuu,2)=Goodness;
XdFINAL(uuu,2)=Goodness;
CTFINAL(uuu,2)=Goodness;
HFINAL(uuu,2)=Goodness;
TpFINAL(uuu,2)=Goodness;
FeFINAL(uuu,2)=Goodness;


subaxis(1,3,1,'SpacingVert',0.04,'SpacingHoriz',0.04)
p2=semilogy(CREEf./norm,'-','Color',[0.5 0.5 0.5]);
p2.Color(4)=Goodness*2;



end

HF=randsample(HFINAL(:,1),3000,true,HFINAL(:,2)./(sum(HFINAL(:,2))));
MeanHF=mean(HF);
StdHF=2*std(HF);

Pen=randsample(Penriched(:,1),3000,true,Penriched(:,2)./(sum(Penriched(:,2))));
MeanPen=mean(Pen);
StdPen=2*std(Pen);

Ur2=randsample(UrinFINAL(:,1),3000,true,UrinFINAL(:,2)./(sum(UrinFINAL(:,2))));
MeanUr=mean(Ur2);
StdUr=2*std(Ur2);

FeF2=randsample(FeFINAL(:,1),3000,true,FeFINAL(:,2)./(sum(FeFINAL(:,2))));
MeanFe=mean(FeF2);
StdFe=2*std(FeF2);

TpF=randsample(TpFINAL(:,1),3000,true,TpFINAL(:,2)./(sum(TpFINAL(:,2))));
MeanTp=mean(TpF);
StdTp=2*std(TpF);

XdF=randsample(XdFINAL(:,1),3000,true,XdFINAL(:,2)./(sum(XdFINAL(:,2))));
MeanXd=mean(XdF);
StdXd=2*std(XdF);

subaxis(1,3,1,'SpacingVert',0.04,'SpacingHoriz',0.04)
errorbar(ResEa./norm,(ResEb)./norm,'->k','MarkerFaceColor','red','MarkerSize',8)

subaxis(1,3,2,'SpacingVert',0.04,'SpacingHoriz',0.04)
histogram(Pen,'Normalization','pdf','NumBins',18,'FaceColor','red')
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'YAxisLocation','left')
xlabel('Proportion of enriched melt','FontSize',16)
ylabel('Probability density','FontSize',16)
grid off
box on

subaxis(1,3,3,'SpacingVert',0.04,'SpacingHoriz',0.04)
histogram(FeF2,'Normalization','pdf','NumBins',18,'FaceColor','red')
hold on
set(gca,'LineWidth',1,'TickLength',[0.01 0.01],'FontName', 'Times New Roman','FontSize',12,'YAxisLocation','left')
xlabel('\delta^{56}Fe','FontSize',16)
ylabel('Probability density','FontSize',16)
grid off
box on
ylim([0 120])

% sort all variables to find best fits
XdFINAL=sortrows(XdFINAL,2,'descend');
Xd=XdFINAL(1,1);
UrinFINAL=sortrows(UrinFINAL,2,'descend');
Ur=UrinFINAL(1,1)-1;
UchangeFINAL=sortrows(UchangeFINAL,2,'descend');
Uchange=UchangeFINAL(1,1);
Dchange=2;
dtopFINAL=sortrows(dtopFINAL,2,'descend');
dtop=dtopFINAL(1,1);
FfactorFINAL=sortrows(FfactorFINAL,2,'descend');
Ffactor=FfactorFINAL(1,1);
TpFINAL=sortrows(TpFINAL,2,'descend');
Tp=TpFINAL(1,1);
FracPyx=1-Xd;

% rerun Melt=-PX for best fit
[Pressure,F_Per,F_Pyx,T,T_solidadiabat,TsolPer,TsolPyx]=MeltPX(Tp,FracPyx,FracCpx,SiO2,TiO2,Al2O3,Cr2O3,FeO,MnO,MgO,CaO,Na2O,K2O);
[Ctfinal,Pdd,CREEf,CREE,CfH2O,Hflux,Urin,PPP,Cmix,CH2OGarnet,Feiso]=MeltPXtraceDirichletFe(Pressure,F_Pyx,F_Per,T,Dsource,Esource,dtop,Ur,Uchange,Dchange,Xd,Ffactor,m,N,dFePdstart,dFePystart);
subaxis(1,3,1,'SpacingVert',0.04,'SpacingHoriz',0.04)
semilogy(CREEf./norm,'-o','Color',[0 0 0],'LineWidth',1,'MarkerFaceColor',[0.6 0.6 0.6]);
