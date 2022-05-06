%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- Coupled PDE Model ----------------------
%------------------------ B. Griffiths -------------------------
%--------------------- University of Oxford --------------------
%------------------------ 08 April 2021 ------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc;
close all;
load('dx.mat');
load('dy.mat');
load('dz.mat');
load('Intensityofint.mat');
%load('Aberationscaling.mat')
load('Defocusphase.mat')
load('logPhotoionisationlookup2.mat')
kpts2=importdata('c.out.kpts.data');
bands2=importdata('c.out.bands.data');
load('DOS.mat');
Intensityofint=Intensityofint./max(max(Intensityofint));
e=1.6e-19;
dE=0.02*e;
Eg=5.41*e:dE:(24.57*e);
Eg2=5.33*e:dE:5.41*e;%(Electron+|Hole|-Binding)
Eg3=10.648*e:dE:10.66*e;%(Exciton+Exciton-Binding)
Eg4=8.908*e:dE:10.648*e;
Ex=0.08.*e;
bx=0.012.*e;
DP=1.74.*e;
lambda = 790e-9; %wavelength of laser used
n2=2.40; % refractive index of sample
c= 3e8; %speed of light
omega=(2.*pi.*c)./lambda; %angular frequency of laser
m_e= 9.11e-31;%mass of electron
m_eeff=0.48.*m_e; %effective electron mass
epsilon0=8.85e-12; %permittivity of free space
h= 6.62e-34; %plancks constant
hbar=h./(2.*pi);
Ephoton=((h.*c)./lambda);
na= 1.77e29;%density of atomic nuclei per m^3
Massatom= 12.*1.66054e-27;%atomic mass
xofint=11; %odd number between 1 and W
zofint=159; %odd number between 1 and Length
Defocusphaseofint=zeros(1,xofint);
W=1650;
for i=1:xofint
Defocusphaseofint(i)=Defocusphase((W/2+1)+i,(W/2+1)+i);
end
clear Defocusphase
D=-(zofint/2)*dz:dz:(zofint/2)*dz; 

%%
dt= 1e-16; %time per step in carriers
NA=1.4;
PD=300;
Timeapprox=1;
PulseE=10.^(-2:0.03:1.2);
Int=PulseE.*(((6.061e17.*4./pi)./Timeapprox).*(NA./1.4).^4.*(30.*Timeapprox./PD));%.*Aberationscaling;
pupilRad = 4.20e-3; %pupil radius in metre (equal to NA * focal length)
NA=1.4;
Depth=20e-6;
t= 0:dt:(15000*dt);
t2=0:dt:(numel(t)+round((omega.*n2.^2.*D(1))./(c.^2.*Defocusphaseofint(xofint).*dt)) -3)*dt;
sigma=((1810./Timeapprox).*(PD./300)).*dt; %Standard Deviation of temporal component of pulse check figure with Patrick %30fs
mu=((10000).*dt); %Mean poition of temprol component of pulse%10000%-4661
t=t';
E_t=(1/(2*pi*(sigma.^2)))*(exp(-(((t-mu).^2)/(2*(sigma.^2))))).*cos(omega.*t); %gaussian %
E_t=E_t.*(1./(1/(2*pi*(sigma.^2))));%normalisation to 1
tsamp=t(9950:10050);
Itsamp=E_t(9950:10050).^2;
figure(1)
hold off
plot(t.*Timeapprox,E_t.^2);
xlabel('Time (ps)')
ylabel('Normalised Intensity')
set(gca,'fontsize', 16);
grid on
xticklabels({'0','0.5','1','1.5','2'})
axes('Position',[.7 .7 .2 .2])
box on
plot(tsamp,Itsamp)
xlim([9950e-16 10050e-16])
xticklabels({0.995 1 1.005})

%%
Paramatervar=PulseE;

Intensityofinttot=zeros(xofint,zofint,numel(Paramatervar));
for i=1:numel(Paramatervar)
for k=1:zofint
for s=1:xofint
Intensityofinttot(s,k,i)=(Intensityofint(s,k)).*Int(i);
end
end
end

Efull=squeeze(zeros(xofint,zofint,numel(t2),numel(Paramatervar)));
for i=-round((omega.*n2.^2.*D(1))./(c.^2.*Defocusphaseofint(xofint).*dt))+1:numel(t2)-1
for j=round(159/2 -(zofint/2)+1):round(159/2 +(zofint/2))
for s=1:xofint
for q=1:numel(PulseE)
Efull(s,j-round(159/2 -(zofint/2)),i,q)=E_t(i-round((omega.*n2.^2.*D(j))./(c.^2.*Defocusphaseofint(s).*dt))).*sqrt((Intensityofinttot(s,j,q).*2)./(c.*n2.*epsilon0));
end
end
end
end
clear Intensityofint
clear Intensityofinttot
clear Int

Ifull=Efull.^2.*c.*n2.*epsilon0.*0.5;

%%
b1=2.*pi./(0.357e-9);
b2=2.*pi./(0.357e-9);
b3=2.*pi./(0.357e-9);
k=(kpts2(:,2).*b1) + (kpts2(:,3).*b2) + (kpts2(:,4).*b3);

bands=bands2-9.938-0.6493;
clear bands2
for j=129:262
for i=1:34
bands(j,i)=bands(j,i)+0.966;
end
end
figure(3) %sort out number
plot(bands')
xticks([0 10 17 21 25 33])
xticklabels({'W','G','X','W','L','G'})
xlabel('K Space')
ylabel('Energy (eV)')
title('Band Structure')

Energy=-25:0.1:25;

DOS=zeros(1,numel(Energy));
kfun=zeros(size(bands,2),numel(Energy));
for x=2:numel(Energy)
for i=1:262
for j=1:34
    if Energy(x-1)<=bands(i,j)
    if bands(i,j)<=Energy(x)
DOS(x)=DOS(x)+1;
kfun(j,x)=kfun(j,x)+1;
    end
    end
end
end
end
DOS=smooth(DOS);
DOS(252:304)=0;
for j=445:numel(DOS)
DOS(j)=DOS(421+444-j)+8;
end
DOS=DOS.*((4.*na)./(sum(DOS(1:250))));

figure(4)
plot(Energy,DOS.*1e-27,'Linewidth',2)
xlabel('Energy (eV)')
ylabel('Density of states (10^{21} cm^{-3})')
set(gca,'fontsize', 16);
grid on
xlim([-25 25])

save('DOS.mat','DOS')

Rhocritcond=zeros(1,numel(Eg));
for j=305:numel(DOS)-1
for i=1+5*(j-305):5*(j-305)+5
Rhocritcond(i)=DOS(j);
end
end
Rhocritcond=Rhocritcond(1:numel(Eg));

Rhocritval=zeros(1,numel(Eg));
for j=1:193
for i=1+5*(j-1):5*(j-1)+5
Rhocritval(i)=DOS(253-j);
end
end
Rhocritval=Rhocritval(1:numel(Eg));


Tau_elecphonon=1./(Rhocritcond.*9.9305e-14);
Tau_holephonon=1./(Rhocritval.*9.9305e-14);
Phononenergy=0.15.*e;
Tau_collph=((1./Tau_elecphonon)+(1./Tau_holephonon))^-1;
Tau_phe=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
Tau_phh=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
Tau_ph=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
dRhodtE=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
dRhodtH=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));

A=1000;
RhoE=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
for j=1:numel(Eg)
RhoE(:,:,A,j,:)=2.807e26./(exp(((Eg(j))-(4.4008.*e))./(1.38e-23.*300))+1);
end
RhoH=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
for j=1:numel(Eg)
RhoH(:,:,A,j,:)=2.807e26./(exp(((Eg(j))-(1.*e))./(1.38e-23.*300))+1);
end
n=squeeze(zeros(xofint,zofint,numel(t2),numel(Paramatervar)));
p=squeeze(zeros(xofint,zofint,numel(t2),numel(Paramatervar)));
n(:,:,A,:)=sum(RhoE(:,:,A,:,:),4);
p(:,:,A,:)=sum(RhoH(:,:,A,:,:),4);

plasmafreq=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
realepsilon=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
imepsilon=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
n3=squeeze(zeros(xofint,zofint,numel(t2),numel(Paramatervar)));
Kappa=squeeze(zeros(xofint,zofint,numel(t2),numel(Paramatervar)));
DeltaTdt=squeeze(zeros(xofint,zofint,numel(t2),numel(Paramatervar)));

Tau_coll=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
Taucoll_av=squeeze(zeros(xofint,zofint,numel(t2),numel(Paramatervar)));
Tau_colleh=squeeze(zeros(xofint,zofint,numel(t2),numel(Paramatervar)));
mobeh=squeeze(zeros(xofint,zofint,numel(t2),numel(Paramatervar)));
Aval=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
Auger=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
dEavdt1=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
dEavdt3=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
Elecphononphoton=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));

averageEE=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
averageEnE=ones(xofint,zofint,numel(t2),numel(Paramatervar));
averageEH=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
averageEnH=ones(xofint,zofint,numel(t2),numel(Paramatervar));

Noption=zeros(1,numel(Eg));
Microstates=zeros(1,numel(Eg));
carrier=zeros(numel(Eg),numel(Eg));
for j=round((5.41.*e)./dE):numel(Eg)
Noption(j)=j-round((5.41.*e)./dE)+1;
Microstates(j)= Microstates(j-1)+(Noption(j));
for x=1:Noption(j)
carrier(j,j-round((5.41.*e)./dE)+2-x)=sum(1./(1+(1:x)));
end
end
carrier(carrier<0)=0;
mult=sum(carrier,2);
carrier=carrier./mult;

AvallossE=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
AvalgainE=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
AvalgainEH=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
AvallossH=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
AvalgainH=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
AvalgainHE=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
AugerlossE=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
AugergainE=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
AugerlossEH=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
AugerlossH=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
AugergainH=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
AugerlossHE=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));

Photoionisation=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
PhononPhotonPropE=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
PhononPhotonPropH=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
PhononPhotonJumpsE=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
PhononPhotonProportionE=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
PhononPhotonJumps2E=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
PhononPhotonProportion2E=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
PhononPhotonJumpsH=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
PhononPhotonProportionH=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
PhononPhotonJumps2H=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
PhononPhotonProportion2H=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));

Temp=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
Temp(:,:,A,:,:)=300;
ETemp=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
ETemp(:,:,A,:,:)=300;
HTemp=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
HTemp(:,:,A,:,:)=300;

Excgain=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg2),numel(Paramatervar)));
Exc=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg2),numel(Paramatervar)));
BiExc=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg2),numel(Paramatervar)));
BiExcgain=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg3),numel(Paramatervar)));
BiExctot=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg3),numel(Paramatervar)));
Exctot=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg2),numel(Paramatervar)));
Dissoc1=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg2),numel(Paramatervar)));
Dissoc2=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg3),numel(Paramatervar)));
Dissoc3=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg4),numel(Paramatervar)));
STbE=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg4),numel(Paramatervar)));
STbEtot=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg4),numel(Paramatervar)));
dFrenk=squeeze(zeros(xofint,zofint,numel(t2),numel(Paramatervar)));
Frenk=squeeze(zeros(xofint,zofint,numel(t2),numel(Paramatervar)));

FD_E=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
FD_H=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));


Ereflect=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
Atten=squeeze(ones(xofint,zofint,numel(t2),numel(Paramatervar)));
PropR=squeeze(zeros(xofint,zofint,numel(t2),numel(Eg),numel(Paramatervar)));
Rsq=zeros(xofint,xofint,zofint,zofint);

ii=2:zofint;
iii=1:xofint;
for k=1:zofint
for s=1:xofint
Rsq(s,iii,k,ii)=((dx.^2)./((dx.*(s-iii)').^2 + (dz.*(k-ii)).^2));
end
end
Rsq(Rsq==Inf)=1;

EN=Eg(1:880)./e;

%%

q=1:numel(Paramatervar);
for i=A:numel(t2) -1
k=1:zofint;
s=1:xofint;
j=1:numel(Eg);
n(:,:,i,:)=sum(RhoE(:,:,i,:,:),4);
p(:,:,i,:)=sum(RhoH(:,:,i,:,:),4);
plasmafreq(s,k,i,q)=sqrt(((n(s,k,i,q)+ p(s,k,i,q)).*e.^2)./(m_eeff.*epsilon0));
mobeh(s,k,i,q)=1e-4.*((4.1e15.*Temp(s,k,i,q).^(3/2))./(((n(s,k,i,q)).*1e-4).*log(1+(5e5.*Temp(s,k,i,q).^2.*((n(s,k,i,q)).*1e-4).^(-(2/3))))));
Tau_colleh(s,k,i,q)=mobeh(s,k,i,q).*m_eeff./(e); %Thermalisation
Tau_coll(s,k,i,j,q)=((1./Tau_collph(j))+(1./(Tau_colleh(s,k,i,q)))).^-1;

Aval(s,k,i,j,q)=((Timeapprox.*dt)./(Tau_colleh(s,k,i,q))).*((2e26-(p(s,k,i,q)./numel(Eg)))./(4.*na)).*(((4.*na)-(n(s,k,i,q)))./(4.*na)).*mult(j)';
Auger(s,k,i,j,q)=(Timeapprox.*dt)./(Tau_colleh(s,k,i,q).*(4.*na).^2).*((n(s,k,i,q).*(p(s,k,i,q)./numel(Eg)))).*mult(j)';
Aval(Aval>1)=1;
Aval(Aval<0)=0;
Auger(Auger>1)=1;
Auger(Auger<0)=0;

k=1:zofint;
for j=1:numel(Eg)
averageEE(s,k,i,j,q)=(j.*RhoE(s,k,i,j,q));
averageEH(s,k,i,j,q)=(j.*RhoH(s,k,i,j,q));
end
averageEnE(s,k,i,q)=round(sum(averageEE(s,k,i,:,q),4)./(n(s,k,i,q)));
averageEnE(averageEnE==0)=1;
averageEnH(s,k,i,q)=round(sum(averageEH(s,k,i,:,q),4)./(p(s,k,i,q)));
averageEnH(averageEnH==0)=1;
for j=1:numel(Eg)
Taucoll_av(s,k,i,q)=Taucoll_av(s,k,i,q)+ (RhoE(s,k,i,j,q)./Tau_coll_e(s,k,i,j,q)) + (RhoH(s,k,i,j,q)./Tau_coll_h(s,k,i,j,q));   
end
Taucoll_av(s,k,i,q)=((1./(n(s,k,i,q)+ p(s,k,i,q))).*Taucoll_av(s,k,i,q)).^-1;
j=1:numel(Eg);

realepsilon(s,k,i,q)=1-((plasmafreq(s,k,i,q).^2.*Taucoll_av(s,k,i,q).^2)./(1+(omega.^2.*Taucoll_av(s,k,i,q).^2)));
imepsilon(s,k,i,q)=((plasmafreq(s,k,i,q).^2.*Taucoll_av(s,k,i,q))./(omega.*(1+(omega.^2.*Taucoll_av(s,k,i,q).^2))));
n3(s,k,i,q)=(sqrt(((sqrt(realepsilon(s,k,i,q).^2+imepsilon(s,k,i,q).^2))+realepsilon(s,k,i,q))./2));
Kappa(s,k,i,q)=(sqrt(((sqrt(realepsilon(s,k,i,q).^2+imepsilon(s,k,i,q).^2))-realepsilon(s,k,i,q))./2));
for s=1:xofint-1
for k=1:zofint
for ii=1:k-2
trip=round(sum(((omega.*n2.^2.*D(ii:k))./(c.^2.*Defocusphaseofint(s).*dt))-((omega.*n2.^2.*D(ii:k))./(c.^2.*Defocusphaseofint(s+1).*dt)))./(((omega.*n2.^2.*D(2))./(c.^2.*Defocusphaseofint(s).*dt))-((omega.*n2.^2.*D(1))./(c.^2.*Defocusphaseofint(s+1).*dt))));
trip(trip>(31-s))=31-s;
trip(s+trip<1)=1-s;
Atten(s,k,i,q)=Atten(s,k,i,q).*exp(-dz.*Kappa(s+trip,ii,i+round((omega.*n2.^2.*(D(ii)-D(k)))./(c.^2.*Defocusphaseofint(s).*dt)),q).*Defocusphaseofint(s)./(n2));
end
end
end

ii=2:zofint;
iii=1:xofint;
PropR(iii,ii,i,q)=(abs((n2.*n3(iii,ii,i,q)-n2.*n3(iii,ii-1,i,q))./(n2.*n3(iii,ii,i,q)+n2.*n3(iii,ii-1,i,q))).^2);
k=1:zofint;
for s=1:xofint
for iii=1:xofint
for ii=2:zofint
for q=1:numel(PulseE)
Ereflect(s,k,i,q)=Ereflect(s,k,i,q)+squeeze(PropR(iii,ii,i,q).*(Efull(iii,ii,i-round(sqrt((dx.*(s-iii)).^2 + (dz.*(ii-k)).^2).*(n2)./(c.*dt)),q).*Rsq(s,iii,k,ii)))';    
end
end
end
end

for s=1:xofint
for k=1:zofint
for q=1:numel(PulseE)
Efull(s,k,i,q)=(Efull(s,k,i-(round((omega.*(n2.*n3(s,k,i,j,q)).^2.*dz)./(c.^2.*Defocusphaseofint(s).*dt)) - round((omega.*(n2.*1).^2.*dz)./(c.^2.*Defocusphaseofint(s).*dt))),q)).*Atten(s,k,i,q) - Ereflect(s,k,i,q);
Ifull(s,k,i,q)=Efull(s,k,i,q).^2.*(c.*n2.*epsilon0).*0.5;
end
end
end
q=1:numel(PulseE);

for k=1:zofint
for s=1:xofint
for j=1:numel(Eg)
Photoionisation(s,k,i,j,q)=(10.^logPhotoionisationlookup2(Efull(s,k,i,q)));%./n3(s,k,i,j,q)
end
end
end
k=1:zofint;
s=1:xofint;
j=1:numel(Eg);
%----------------
dEavdt1(s,k,i,q)=((e.^2.*Ifull(s,k,i,q))./(m_eeff.*c.*n2.*n3(s,k,i,q).*epsilon0));
dEavdt3(s,k,i,j,q)=(((omega).^2 .* (Tau_coll(s,k,i,j,q)).^2)+1);
Elecphononphoton(s,k,i,j,q)=(Tau_coll(s,k,i,j,q).*dEavdt1(s,k,i,j,q)./dEavdt3(s,k,i,j,q));

for j=1:numel(Eg)
PhononPhotonPropE(s,k,i,j,q)=((Elecphononphoton(s,k,i,q).*dt)./(Ephoton-Phononenergy));
end

for k=1:zofint
a=squeeze(log10(RhoE(s,k,i,1:880,q)));
[fitresult] = ETempfit(EN, a);
coeffvals= coeffvalues(fitresult);
FD_E(s,k,i,1:880,q)=10.^(coeffvals(1).*EN(:) +coeffvals(2));
ETemp(s,k,i,q)=(-(Eg(300)-Eg(30)))./(1.38e-23.*log((10.^(coeffvals(1).*EN(300) +coeffvals(2)))./(10.^(coeffvals(1).*EN(30) +coeffvals(2)))));
a=squeeze(log10(RhoH(s,k,i,1:880,q)));
[fitresult] = ETempfit(EN, a);
coeffvals= coeffvalues(fitresult);
FD_H(s,k,i,1:880,q)=10.^(coeffvals(1).*EN(:) +coeffvals(2));
HTemp(s,k,i,q)=(-(Eg(300)-Eg(30)))./(1.38e-23.*log((10.^(coeffvals(1).*EN(300) +coeffvals(2)))./(10.^(coeffvals(1).*EN(30) +coeffvals(2)))));
end
k=1:zofint;
FD_E(s,k,i,1:880,q)=n(s,k,i,q).*(FD_E(s,k,i,1:880,q)./sum(FD_E(s,k,i,1:880,q),4));
FD_H(s,k,i,1:880,q)=p(s,k,i,q).*(FD_H(s,k,i,1:880,q)./sum(FD_H(s,k,i,1:880,q),4));


for j=round((5.41.*e)./dE):numel(Eg)
for x=1:Noption(j)
a=(RhoE(s,k,i,x,q) + (RhoE(s,k,i,j,q).*Aval(s,k,i,j,q).*(2.*carrier(j,x))))<Rhocritcond(j);
AvallossE(s,k,i,j,q)=AvallossE(s,k,i,j,q)+(a.*Aval(s,k,i,j,q).*RhoE(s,k,i,j,q).*carrier(j,x));
AvalgainE(s,k,i,x,q)=AvalgainE(s,k,i,x,q)+(Aval(s,k,i,j,q).*RhoE(s,k,i,j,q).*(2.*carrier(j,x)));
AvalgainEH(s,k,i,x,q)=AvalgainEH(s,k,i,x,q)+(Aval(s,k,i,j,q).*RhoE(s,k,i,j,q).*(carrier(j,x)));
a=(RhoH(s,k,i,x,q) + RhoH(s,k,i,j,q).*Aval(s,k,i,j,q).*(2.*carrier(j,x)))<Rhocritval(x);
AvallossH(s,k,i,j,q)=AvallossH(s,k,i,j,q)+(a.*Aval(s,k,i,j,q).*RhoH(s,k,i,j,q).*carrier(j,x));
AvalgainH(s,k,i,x,q)=AvalgainH(s,k,i,x,q)+(a.*Aval(s,k,i,j,q).*RhoH(s,k,i,j,q).*(2.*carrier(j,x)));
AvalgainHE(s,k,i,x,q)=AvalgainHE(s,k,i,x,q)+(a.*Aval(s,k,i,j,q).*RhoH(s,k,i,j,q).*(carrier(j,x)));
a=(RhoE(s,k,i,j,q) + RhoE(s,k,i,x,q).*Auger(s,k,i,j,q).*carrier(j,x)) <Rhocritcond(j);
AugerlossE(s,k,i,x,q)=AugerlossE(s,k,i,x,q)+(a.*Auger(s,k,i,j,q).*RhoE(s,k,i,x,q).*(2.*carrier(j,x)));
AugerlossEH(s,k,i,x,q)=AugerlossEH(s,k,i,x,q)+(a.*Auger(s,k,i,j,q).*RhoE(s,k,i,x,q).*(carrier(j,x)));
a=(RhoH(s,k,i,j,q) + RhoH(s,k,i,x,q).*Auger(s,k,i,j,q).*carrier(j,x)) <Rhocritval(j);
AugerlossH(s,k,i,x,q)=AugerlossH(s,k,i,x,q)+(a.*Auger(s,k,i,j,q).*RhoH(s,k,i,x,q).*(2.*carrier(j,x)));
AugerlossHE(s,k,i,x,q)=AugerlossHE(s,k,i,x,q)+(a.*Auger(s,k,i,j,q).*RhoH(s,k,i,x,q).*(carrier(j,x)));
end
end
for j=round((5.41.*e)./dE):numel(Eg)
AugergainE(s,k,i,j,q)=0.5.*AugerlossE(k,i,j - round((5.41.*e)./dE) +1);
AugergainH(s,k,i,j,q)=0.5.*AugerlossH(k,i,j - round((5.41.*e)./dE) +1);
PhononPhotonPropE(:,:,i,j,q)=((Elecphononphoton(:,:,i,j,q).*dt)./(Ephoton-Phononenergy));
PhononPhotonPropH(:,:,i,j,q)=(Holephononphoton(:,:,i,j,q).*dt)./(Ephoton-Phononenergy);
end
PhononPhotonJumpsE(:,:,i,:,:)=ceil(PhononPhotonPropE(:,:,i,:,:));
PhononPhotonJumpsE(PhononPhotonJumpsE<1)=1;
PhononPhotonProportionE(:,:,i,:,:)=PhononPhotonPropE(:,:,i,:,:)-floor(PhononPhotonPropE(:,:,i,:,:));
PhononPhotonJumps2E(:,:,i,:,:)=PhononPhotonJumpsE(:,:,i,:,:)-1;
PhononPhotonProportion2E(:,:,i,:,:)=1-PhononPhotonProportionE(:,:,i,:,:);
PhononPhotonJumpsH(:,:,i,:,:)=ceil(PhononPhotonPropH(:,:,i,:,:));
PhononPhotonJumpsH(PhononPhotonJumpsH<1)=1;
PhononPhotonProportionH(:,:,i,:,:)=PhononPhotonPropH(:,:,i,:,:)-floor(PhononPhotonPropH(:,:,i,:,:));
PhononPhotonJumps2H(:,:,i,:,:)=PhononPhotonJumpsH(:,:,i,:,:)-1;
PhononPhotonProportion2H(:,:,i,:,:)=1-PhononPhotonProportionH(:,:,i,:,:);

%Carrier_phonon_photon
for j=1:numel(Eg)-round(Ephoton./dE)
if j+(PhononPhotonJumpsE(s,k,i,j,q)*round((Ephoton-Phononenergy)./dE))<numel(Eg)
RhoE(s,k,i+1,j+(PhononPhotonJumpsE(s,k,i,j,q).*round((Ephoton-Phononenergy)./dE)),q)=RhoE(s,k,i+1,j+(PhononPhotonJumpsE(s,k,i,j,q).*round((Ephoton-Phononenergy)./dE)),q) + (PhononPhotonProportionE(s,k,i,j,q).*RhoE(s,k,i,j,q));
RhoE(s,k,i+1,j+(PhononPhotonJumps2E(s,k,i,j,q).*round((Ephoton-Phononenergy)./dE)),q)=RhoE(s,k,i+1,j+(PhononPhotonJumps2E(s,k,i,j,q).*round((Ephoton-Phononenergy)./dE)),q) + (PhononPhotonProportion2E(s,k,i,j,q).*RhoE(s,k,i,j,q));
end
if j+(PhononPhotonJumpsH(s,k,i,j,q)*round((Ephoton-PhononenergyH(j))./dE))<numel(Eg)
RhoH(s,k,i+1,j+(PhononPhotonJumpsH(s,k,i,j,q).*round((Ephoton-PhononenergyH(j))./dE)),q)=RhoH(s,k,i+1,j+(PhononPhotonJumpsH(s,k,i,j,q).*round((Ephoton-PhononenergyH(j))./dE)),q) + (PhononPhotonProportionH(s,k,i,j,q).*RhoH(s,k,i,j,q));
RhoH(s,k,i+1,j+(PhononPhotonJumps2H(s,k,i,j,q).*round((Ephoton-PhononenergyH(j))./dE)),q)=RhoH(s,k,i+1,j+(PhononPhotonJumps2H(s,k,i,j,q).*round((Ephoton-PhononenergyH(j))./dE)),q) + (PhononPhotonProportion2H(s,k,i,j,q).*RhoH(s,k,i,j,q));
end
end

Tau_phe(s,k,i,j,q)=Tau_elecphonon(j).*((ETemp(s,k,i,q)-Temp(s,k,i,q))./(ETemp(s,k,i,q)+Temp(s,k,i,q)));
Tau_phh(s,k,i,j,q)=Tau_holephonon(j).*((HTemp(s,k,i,q)-Temp(s,k,i,q))./(HTemp(s,k,i,q)+Temp(s,k,i,q)));
Tau_ph(s,k,i,j,q)=Tau_collph(j).*((ETemp(s,k,i,q)-Temp(s,k,i,q))./(ETemp(s,k,i,q)+Temp(s,k,i,q)));

for j=8:numel(Eg)
RhoE(s,k,i+1,j,q)=RhoE(s,k,i+1,j,q) - ((dt./abs(Tau_phe(s,k,i,j,q))).*RhoE(s,k,i,j,q));
RhoH(s,k,i+1,j,q)=RhoH(s,k,i+1,j,q) - ((dt./abs(Tau_phh(s,k,i,j,q))).*RhoH(s,k,i,j,q));
RhoE(s,k,i+1,j-sign((dt./Tau_phe(s,k,i,j,q)).*(round(Phononenergy./dE)),q))=RhoE(s,k,i+1,j-sign((dt./Tau_phe(s,k,i,j,q)).*(round(Phononenergy./dE)),q)) + ((dt./abs(Tau_phe(s,k,i,j,q))).*RhoE(s,k,i,j,q));
RhoH(s,k,i+1,j-sign((dt./Tau_phh(s,k,i,j,q)).*(round(Phononenergy./dE)),q))=RhoH(s,k,i+1,j-sign((dt./Tau_phh(s,k,i,j,q)).*(round(Phononenergy./dE)),q)) + ((dt./abs(Tau_phh(s,k,i,j,q))).*RhoH(s,k,i,j,q));
end

a=Tau_colleh(:,i);
a(a<1e-15)=1.1e-15;
dRhodtE(s,k,i,1:880,q)=(FD_E(s,k,i,1:880,q)-RhoE(s,k,i,1:880,q))./(a);
dRhodtH(s,k,i,1:880,q)=(FD_H(s,k,i,1:880,q)-RhoH(s,k,i,1:880,q))./(a);
for j=1:numel(Eg)
RhoE(s,k,i+1,j,q)=RhoE(s,k,i+1,j,q)+(dRhodtE(s,k,i,j,q).*dt.*Timeapprox)-(AvallossE(s,k,i,j,q))+(AvalgainE(s,k,i,j,q)) +(AvalgainHE(s,k,i,j,q))-(AugerlossE(s,k,i,j,q))+(AugergainE(s,k,i,j,q))-(AugerlossHE(s,k,i,j,q))+ Photoionisation(s,k,i,j,q)*dt - (RhoE(s,k,i+1,j,q).*dt.*Timeapprox./carrierradtau);% +((gradDiffE(:,i).*(RhoE(s,k,i,j,q)./n(:,i))).*dt.*Timeapprox);
RhoH(s,k,i+1,j,q)=RhoH(s,k,i+1,j,q)+(dRhodtH(s,k,i,j,q).*dt.*Timeapprox)+(AvalgainEH(s,k,i,j,q)) -(AvallossH(s,k,i,j,q))+(AvalgainH(s,k,i,j,q))-(AugerlossH(s,k,i,j,q))+(AugergainH(s,k,i,j,q))-(AugerlossEH(s,k,i,j,q)) - (RhoH(s,k,i+1,j,q).*dt.*Timeapprox./carrierradtau);%+ ((gradDiffH(:,i).*(RhoH(s,k,i,j,q)./p(:,i))).*dt.*Timeapprox);
end
RhoH(:,:,i+1,1,:)=RhoH(:,:,i+1,1,:)+sum(Photoionisation(:,:,i,:,:),4)*dt;
n(:,:,i+1,:)=sum(RhoE(:,:,i+1,:,:),4);
p(:,:,i+1,:)=sum(RhoH(:,:,i+1,:,:),4);

for j=1:numel(Eg)
DeltaTdt(:,:,i+1,:)=DeltaTdt(:,:,i+1,:)+((1.38e-23.*Phononenergy.*(RhoE(s,k,i+1,j,q)))./(na.*Tau_phe(s,k,i,j,q)))+((1.38e-23.*Phononenergy.*(RhoH(s,k,i+1,j,q)))./(na.*Tau_phh(s,k,i,j,q)));
end
Temp(:,:,i+1,:) = Temp(:,:,i,:) + DeltaTdt(:,:,i+1,:).*dt.*Timeapprox;
k=1:zofint;

%Excitons

for j=1:numel(Eg2)
x=1:numel(Eg2);
y=round(((Eg(x)+Eg(j)-(Eg(1))-Ex) - Eg2(1))./(dE))+1;
x(y>numel(Eg2))=[];
y(y>numel(Eg2))=[];
Excgain(s,k,i+1,y,q)=(((RhoE(s,k,i+1,j,q).*RhoH(s,k,i+1,x,q)).*Timeapprox.*dt)./(na.*Tau_colleh(s,k,i,q)));
RhoE(s,k,i+1,j,q)=RhoE(s,k,i+1,j,q)-(sum(Excgain(s,k,i,y,q),4)./2);
RhoH(s,k,i+1,j,q)=RhoH(s,k,i+1,j,q)-(sum(Excgain(s,k,i,y,q),4)./2);
Exc(s,k,i+1,j,q)=Exc(s,k,i,y,q)+Excgain(s,k,i,y,q);
end
Exctot(s,k,i+1,:,q)=Exctot(s,k,i,:,q)+Exc(s,k,i,:,q) - (Exctot(s,k,i,:,q).*dt.*Timeapprox./Excradtau);
for j=1:numel(Eg2)
Dissoc1(s,k,i,j,q)=PhononPhotonPropE(s,k,i,j,q) +(exp(-(((Eg2(1)+0.08.*e-Eg2(j)))./(1.38e-23.*ElecTemp(s,k,i,q))))).*dt.*Timeapprox./Tau_colleh(s,k,i,q);
end
Dissoc1(Dissoc1>1)=1;
for j=1:numel(Eg2)
RhoE(s,k,i+1,round(((Eg2(j)./e)- Eg(1)./e)./(dE./e)./2)+3,q)=RhoE(s,k,i+1,round(((Eg2(j)./e) -Eg(1)./e)./(dE./e)./2)+3,q) + (Dissoc1(s,k,i,j,q).*Exctot(s,k,i+1,j,q)./2);
RhoH(s,k,i+1,round(((Eg2(j)./e)- Eg(1)./e)./(dE./e)./2)+3,q)=RhoH(s,k,i+1,round(((Eg2(j)./e) -Eg(1)./e)./(dE./e)./2)+3,q) + (Dissoc1(s,k,i,j,q).*Exctot(s,k,i+1,j,q)./2);
Exctot(k,i+1,j)=Exctot(k,i+1,j)-(Dissoc1(k,i,j).*Exctot(k,i+1,j));
end
Exctot(Exctot<0)=0;

for j=1:numel(Eg2)
x=1:numel(Eg2);
y=round(((Eg2(x)+Eg2(j)-Ex) - Eg3(1))./(dE))+1;
x(y>numel(Eg3))=[];
y(y>numel(Eg3))=[];
BiExcgain(s,k,i+1,y,q)=((Exctot(s,k,i+1,j,q).*Exctot(s,k,i+1,x,q).*Timeapprox.*dt)./(na.*Tau_colleh(s,k,i,q)));
Exctot(s,k,i+1,j,q)=Exctot(s,k,i+1,j,q)-sum(BiExcgain(s,k,i+1,y,q),3);
BiExc(s,k,i,y,q)=BiExc(s,k,i,y,q)+BiExcgain(s,k,i,y,q);
end
BiExctot(s,k,i+1,:,q)=BiExctot(s,k,i,:,q)+BiExc(s,k,i,:,q) - (BiExctot(s,k,i,:,q).*dt.*Timeapprox./BiExcradtau);
for j=1:numel(Eg3)
Dissoc2(s,k,i,j,q)=PhononPhotonPropE(s,k,i,j,q) + (exp(-(((Eg3(1)+0.012.*e-Eg3(j)))./(1.38e-23.*ETemp(s,k,i,q))))).*dt.*Timeapprox./Tau_colleh(s,k,i,q);
end
for j=1:numel(Eg3)
Exctot(s,k,i+1,round(((Eg3(j)./e)./2 -Eg2(1)./e)./(dE./e))+3,q)=Exctot(s,k,i+1,round(((Eg3(j)./e)./2 -Eg2(1)./e)./(dE./e))+3,q) + (Dissoc2(s,k,i,j,q).*BiExctot(s,k,i+1,j,q));
BiExctot(s,k,i+1,j,q)=BiExctot(s,k,i+1,j,q)-(Dissoc2(s,k,i,j,q).*BiExctot(s,k,i+1,j,q));
end
BiExctot(BiExctot<0)=0;

for j=1:numel(Eg4)
Dissoc3(s,k,i,j,q)=(exp(-(((Eg4(1)+1.74.*e-Eg4(j)))./(1.38e-23.*ETemp(s,k,i,q))))).*dt.*Timeapprox./Tau_colleh(s,k,i,q);
STbE(s,k,i,j,q)=BiExctot(s,k,i+1,j,q).*dt.*Timeapprox./(Tau_collph(j).*(((dt./(Tau_ph(s,k,i,j,q)))+1)./2));    
end
STbEtot(s,k,i+1,:,q)=STbEtot(s,k,i,:,q)+STbE(s,k,i,:,q)-(STbEtot(s,k,i,:,q).*Dissoc3(s,k,i,:,q));
BiExctot(s,k,i+1,:,q)=BiExctot(s,k,i+1,:,q)- STbE(s,k,i,:,q);

for j=1:numel(Eg4)
dFrenk(s,k,i,q)=STbEtot(s,k,i,j,q).*dt.*Timeapprox.*(((Phononenergy./6.626e-34)).*exp(-((((Eg4(1)+0.47.*e-Eg4(j)))))./(1.38e-23.*Temp(s,k,i,q))));
end
STbEtot(s,k,i+1,q)=STbEtot(s,k,i+1,q) - (dFrenk(s,k,i,j,q))- STbEtot(s,k,i+1,q).*exp(-(1.74.*e)./(1.38e-23.*Temp(s,k,i,j,q))) - STbEtot(s,k,i+1,q).*dt.*Timeapprox./5e-9;

Frenk(:,:,i+1,:)=Frenk(:,:,i,:) + dFrenk(:,:,i,:);
end
save('RhoE.mat','RhoE');
save('STbEtot.mat','STbEtot');
save('Exctot.mat','Exctot');
save('BiExctot.mat','BiExctot');
save('Frenk.mat','Frenk');
save('Efull.mat','Efull');
save('ETemp.mat','ETemp');
save('Temp.mat','Temp');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------- END ------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
