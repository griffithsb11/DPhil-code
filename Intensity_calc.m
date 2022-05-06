close all;
clear, clc;
PulseEnergy=1e-9;

L=825; %length and width of pixels)Y = no of rows
W=825;%X = no. of columns
Length = 1030;%206; %number of z planes included in the calculation of the PSF stack
R = 50; %radius of circular aperture in pixels - this corresponds to the pupil of the objective(:,:,i). It needs to be sufficiently small relative(:,:,i) to L & W to have(:,:,i) reasonable resolution within the PSF
pupilRad = 4.20e-3; %pupil radius in metre (equal to NA * focal length)
lambda = 790e-9; %wavelength of laser used
% f = 200/1000; %focal distance of lens
f = 3e-3; %focal distance of lens (Mag divided by tube lens length of manufacturer)
n1=1.52; %refractive(:,:,i) index of oil (immersion medium of objective(:,:,i) lens)
n2=2.40; % refractive(:,:,i) index of sample
NA = 1.4;
C1 = ((2*pi)/(lambda));


% z= fspecial('gaussian', [5 5], 0.55); % matched to SLM smoothing effects (not always used)
c= 3e8; %speed of light
omega=(2.*pi.*c)./lambda; %angular frequency of laser
e=1.6e-19; %Charge of electron
m_e= 9.11e-31;%mass of electron
m_eeff=0.4.*m_e; %effective(:,:,i) electron mass
m_exciton= m_eeff./2; %effective(:,:,i) exciton mass =1/m_exciton = 1/me + 1/mh
epsilon0=8.85e-12; %permittivity of free space
h= 6.62e-34; %plancks constant
Ephoton= (h.*c)./ lambda; %photon energy
hbar=h./(2.*pi);
na= 1.77e29;%density of atomic nuclei per m^3
Massatom= 12.*1.66054e-27;%atomic mass
xofint=31; %odd number between 1 and W           31
zofint=159; %odd number between 1 and Length         159
%PARAMETERS TO CHANGE



Depth = 20e-6; % depth focused to inside sample -looking at spherical aberration

Depth = double(Depth);

%Here we create the pupil at the back of the objective lens and define a
%top-hat function across it. Our resolution of this defines the resolution
%of of focal pixels.
%PUPILMASK CIRCULAR%
pupilmask = zeros(L,W,'single'); 
m=size(pupilmask,1); %number of rows (= L)%
n=size(pupilmask,2); %number of columns (= W)
    if mod(m,2)==1
        x=(0:m-1)-(m-1)/2;                % odd num
    else
        x=(1:m)-(m+1)/2;                   % eve(:,:,i)n
    end
    
    if mod(n,2)==1
        y=(0:n-1)-(n-1)/2;                 % odd num
    else
        y=(1:n)-(n+1)/2;                  % even
    end

Y=meshgrid(x,y); %[X,Y]
X=Y';Y=X'; %takes the tranpose of the matrixclc and creates a grid of x and y co-ordinates centred at the middle of the pupil mask
% pupilmask(200:300)=1;

X = double(X)*pupilRad/R; %Scaling coordinates in pupil plane
Y=single(Y)*pupilRad/R;
r=(X.^2+Y.^2)/(pupilRad*pupilRad);
%r= pupilmask.*r;
pupilmask(sqrt(X.^2+Y.^2)<(pupilRad))= 1 ; % Top hat intensity distribution in the pupil plane
%pupilmask = exp(-r/0.2); % Gaussian intensity distribution in the pupil plane
pupilmask(sqrt(X.^2+Y.^2)>(pupilRad))= 0;
%pupilmask(sqrt(X.^2+Y.^2)<(19*pupilRad/20))= 0;

r(sqrt(X.^2+Y.^2)>(pupilRad))=0;
pupilmask = single(pupilmask);

Npupil = size((find(pupilmask)),1); 

%calculating dimensions in the focal plane
thetaX = lambda/(2*pupilRad/(R));
thetaY = lambda/(2*pupilRad/(R));
dx = (f*thetaX)/(W*0.5); %physical X distance corresponding to 1 pixel
dy = f*thetaY/(L*0.5);%physical Y distance corresponding to 1 pixel
save('dx.mat','dx')
save('dy.mat','dy')
 

%Define a defocus phase that occurs when we ar not in the plan of focus
Defocusphase =  (C1*NA*sqrt(-(r) +(n2/NA)^2));             %defocus phase, to shift in z in units of metre - Eqn 7 in Jesacher 2010
Defocusphase = pupilmask.*Defocusphase;                                  % Remove phase values outside pupil       
Defocusphase = double(Defocusphase);
DPprime = (sum(sum((Defocusphase(:,:))))/Npupil)*pupilmask;                   % Average defocus in pupil
pupilmask = double(pupilmask);
phi  = -(C1*NA*(sqrt((-1.0*r) +(n2/NA)^2) - sqrt((-1.0*r) + (n1/NA)^2)));         % Phase for Complete RI mismatch compensation
phi = double(phi);
phip = pupilmask.*phi;                                                                      %Only retain phase values in pupil
clear phi;
PhiPprime = (sum(sum((phip(:,:))))/Npupil)*pupilmask;                           %Average SA phase within pupil

Prod1 = (phip - PhiPprime).*(Defocusphase - DPprime);                           %Nominator of Eqn 9 Jesacher 2010 before ave(:,:,i)rage
Prod1Sum = sum(sum((Prod1(:,:))))/Npupil;                                       %Nominator of Eqn 9 Jesacher 2010 
clear Prod1;
Prod2 = (Defocusphase - DPprime).*(Defocusphase - DPprime);                     %deNominator of Eqn 9 Jesacher 2010 before ave(:,:,i)rage                                     
Prod2Sum = sum(sum((Prod2(:,:))))/Npupil;                                       %deNominator of Eqn 9 Jesacher 2010
clear Prod2;
clear PhiPprime;
clear DPprime;

phip = phip - (Prod1Sum/Prod2Sum)*Defocusphase;                                 %Spherical aberration phase minus defocus - only focus distortion        
%phip = X.*pupilmask*pi*1e5;
phip = phip - (sum(sum((phip(:,:))))/Npupil)*pupilmask; % Centring the phase distribution on zero

% Strehl = ((sum(sum(exp(1i*(Depth*phip(:,:)))))) -  (L*W-Npupil))/(Npupil);  %Strehl ratio using the amount of phase aberration
% Strehl = Strehl*conj(Strehl);

%   pupilmask((X)>0.4e-3)=0;  % Slit illumination of the objective pupil 
%   pupilmask((X)<-0.4e-3)=0;
% Npupil = size((find(pupilmask)),1); 
  

%plot(pupilmask(400,301:500), 'DisplayName', 'pupilmask(400,301:500)', 'YDataSource','pupilmask(400,301:500)'); figure(gcf)

clear Rsq
clear Size_pixel_focalplane
clear Test
clear thetaX
clear thetaY
clear x
clear X
clear X4
clear Y
clear y
clear z
clear m
clear n
clear n1
clear r
clear k
clear NA

%pupilmask2 = conv2((pupilmask),z,'same');    %low pass filter that mimics the smoothing of the SLM between adjacent pixels (phase wrapping in particular)


Pupilphase= -(15e-6.*phip); %=zeros(L,W);  %Phase needed to be applied to SLM to account for spherical abberation from depth into sample

save('Defocusphase.mat','Defocusphase');
figure(1)
pcolor(squeeze(Pupilphase(:,:))), shading interp; colorbar
xlim([W/2-65 W/2+65]) % xlim([W/2 W/2+22])
ylim([W/2-65 W/2+65]) % ylim([W/2 W/2+22])
xlabel('Radius')
ylabel('Radius')
title('Spherical Aberation correction applied to SLM')
% savefig('Aberations.fig')

%gaussi = 1.0;%0.5 + 0.2*gauss;
%pupilmask0 = pupilmask.*exp(-r/gaussi); % Gaussian intensity distribution in the pupil plane
D=-20.0e-6; %axial distance of slice furthest from the nominal focal plane. The axial dimensions of the calculation spread from  "Depth + D" to "Depth +(Length*Dstep + D)"
dz = 40.0000e-06./Length; % axial resolution of the calculation
save('dz.mat','dz')
% D = double(D);

Intensity=zeros(L,W,Length);
G2=zeros(L,W,Length);
GX4=zeros(L,W,Length);
    for j = 1:(Length)
    G2(:,:,j) = (abs(pupilmask).*exp(1i*((Depth)*phip(:,:))).*exp(1i*D*Defocusphase(:,:)).*exp(1i*Pupilphase(:,:)));
    GX4(:,:,j)=fftshift(fft2(fftshift(G2(:,:,j)),L,W));
    Intensity(:,:,j) = GX4(:,:,j).*conj(GX4(:,:,j));
    D=D+dz;
    end

% G = (abs(pupilmask).*exp(1i*Depth*phip(:,:));
% X4=fftshift(fft2(fftshift(G),L,W)); 
% Nocorrection = (max(max(X4(:,:))))*conj(max(max(X4(:,:)))); %Intensity given no SLM correction 
G = abs(pupilmask);
X4=fftshift(fft2(fftshift(G),L,W)); 
Noaberration = (max(max(X4(:,:))))*conj(max(max(X4(:,:))));%Intensity given perfect SLM correction
Aberationscaling=max(max(max(Intensity)))/Noaberration;

Aberphase(:,:)=((Depth)*phip(:,:))+Pupilphase(:,:);

clear C1
clear C2
clear D
clear Dstep
clear f
clear G
clear G2
clear GX4
clear X4
clear phip
clear Phip
clear pupilmask
clear Npupil
clear Prod1Sum
clear Prod2Sum
clear Pupilphase
clear pupilRad


Intensity=(abs(Intensity)./max(max(max(Intensity)))).*Aberationscaling;


Scalez_rz=zeros(W,Length);
Scalez_rz(:,1)=20e-6+((Length/2)-1)*dz;
for B=2:Length
Scalez_rz(:,B)=Scalez_rz(:,B-1)-dz;
end

Scaler_rz=zeros(W,Length);
Scaler_rz(1,:)=-W.*dx./2;
for B=2:W
Scaler_rz(B,:)=Scaler_rz(B-1,:)+dx;
end

Iofint=squeeze(Intensity(:,round(W./2),:))';
Iofint=(Iofint./max(max(Iofint))).*Aberationscaling;

figure(2)
pcolor(Scaler_rz',Scalez_rz',squeeze(Iofint(:,:))), shading interp; colorbar
xlabel('Radius (\mum)')
ylabel('Depth (\mum)')
yticks([0, 0.5e-5, 1e-5, 1.5e-5, 2.0e-5, 2.5e-5, 3.0e-5, 3.5e-5,3.99e-05])
xticks([-1.4107e-05, -1e-5, -0.5e-5, 0, 0.5e-5, 1e-5, 1.4073e-05])
yticklabels({40 35 30 25 20 15 10 5 0})
xticklabels({-15 -10 -5 0 5 10 15 20})

% 
% Intensityofint=zeros(xofint,zofint);
% for i=1:xofint
% Intensityofint(i,:)=Intensity((W/2+1)+i,(L/2)+1,(((Length/2)+1)-((zofint/2)-0.5)):(((Length/2)+1)+((zofint/2)-0.5)));
% end
% 
% clear Intensity
% clear I_peakpixel
% clear Power
% clear Pulse_dur
% clear Sum_Intspace
% 
% save('Intensityofint.mat', 'Intensityofint')
% save('Aberphase.mat','Aberphase')
% save('Aberationscaling.mat','Aberationscaling')
% figure(2)
% pcolor(squeeze(Intensityofint(:,:))), shading interp; colorbar
% % Intensityofint=(Intensityofint./max(max(Intensityofint))).*Aberationscaling;

% load('Intensityofint.mat');
% load('Defocusphase.mat')
% load('Aberphase.mat','Aberphase')
% 
% Defocusphaseofint=zeros(1,xofint);
% for i=1:xofint
% Defocusphaseofint(i)=Defocusphase((W/2+1)+i,(W/2+1)+i);
% end
% clear Defocusphase
% D=-(zofint/2)*dz:dz:(zofint/2)*dz;
% Aberphaseofint=zeros(1,xofint);
% for i=1:xofint
% Aberphaseofint(i)=Aberphase((W/2+1)+i,(W/2+1)+i);
% end
% 
% phaseofint=Defocusphaseofint+Aberphaseofint;
% 
% dt= 1e-16; %time per step in carriers
% PD=300;
% Timeapprox=1;
% t= 0:dt:(10000*dt);
% t2=0:dt:(numel(t)+round((omega.*n2.^2.*D(1))./(c.^2.*Defocusphaseofint(xofint).*dt)) -3)*dt;
% sigma=((1810./Timeapprox).*(PD./300)).*dt; %Standard Deviation of temporal component of pulse check figure with Patrick %30fs
% mu=((5000).*dt); %Mean poition of temprol component of pulse%10000%-4661
% t=t';
% E_t=(1/(2*pi*(sigma.^2)))*(exp(-(((t-mu).^2)/(2*(sigma.^2))))).*cos(omega.*t); %gaussian %
% E_t=E_t.*(1./(1/(2*pi*(sigma.^2))));%normalisation to 1
% figure(3)
% hold off
% plot(t.*Timeapprox,E_t.^2);
% xlabel('Time (ps)')
% ylabel('Normalised Intensity')
% set(gca,'fontsize', 16);
% grid on
% 
% Efull=squeeze(zeros(xofint,zofint,numel(t2)));
% for i=-round((omega.*n2.^2.*D(1))./(c.^2.*phaseofint(xofint).*dt))+1:numel(t2)-1
% for j=round(159/2 -(zofint/2)+1):round(159/2 +(zofint/2))
% for s=1:xofint
% Efull(s,j-round(159/2 -(zofint/2)),i)=E_t(i-round((omega.*n2.^2.*D(j))./(c.^2.*phaseofint(s).*dt))).*sqrt((Intensityofint(s,j)));
% end
% end
% end
% 
% Scalez_rz=zeros(xofint,zofint);
% Scalez_rz(:,1)=20e-6-((zofint/2)-1)*dz;
% for B=2:zofint
% Scalez_rz(:,B)=Scalez_rz(:,B-1)+dz;
% end
% 
% Scaler_rz=zeros(xofint,zofint);
% Scaler_rz(1,:)=0;
% for B=2:xofint
% Scaler_rz(B,:)=Scaler_rz(B-1,:)+dx;
% end
% 
% Scalerad=zeros(xofint,xofint);
% Scalerad(1,:)=0;
% for B=2:xofint
% Scalerad(B,:)=Scalerad(B-1,:)+dx;
% end
% 
% figure(4)
% pcolor(Scalez_rz,Scaler_rz,squeeze(Efull(:,:,mu./dt).^2)), shading interp; colorbar
% xlabel('Depth (\mum)')
% ylabel('Radius (nm)')
% % xticks([1 80.5 159])
% xticklabels({17 18 19 20 21 22 23})
% yticklabels({0 50 100 150 200 250 300 350 400 450 500})


% for i=4000:numel(t2)
% figure(2);
% pcolor(squeeze(Efull(:,:,i).^2)), shading interp; colorbar%.*real(exp(1i.*omega.*t(i))).^2))
% caxis([0 1])
% pause(0.05)
% i
% end
% save('Efull.mat','Efull')

