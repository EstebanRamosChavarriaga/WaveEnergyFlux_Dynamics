%% WAVE ENERGY FLUX VARIABILITY ALONG CONTINENTAL SHELVES FROM WAVE CLIMATE DATA-DRIVEN ANALYTICS. RAMOS-CHAVARRIAGA ET AL., 2021. (Esteban Ramos-Chavarriaga 2021)
% This code is free access and all citations on this paper are required.
% Pre-processing stage: run this section before each Part.
% Part 1: Wave Climate variability and statistical trends (Classical description and exploratory analysis)
% Part 2: Coastal wave energy flux availability regimes in the continental shelf (Classification and description)
% Part 3: Climate periodicities and spectral correlation and coherence to the Oceanic Niño Index (Description)
%% Pre-processing: RUN THIS SECTION BEFORE RUNING EITHER PART 1, 2 OR 3.
clear all, close all, clc

% Data acquisition WW3
filebase = 'C:\Users\name\Desktop\pacifico_sur_nort_';

% Data acquisition - ONI
data_ERSST=xlsread('C:\Users\name\Desktop\ERSST');
ONI=data_ERSST(:,5); % ONI data from ERSST v5

% Data acquisition - Bathymetric data https://www.gmrt.org/GMRTMapTool/
[longt, latt, Z] = grdread2('C:\Users\name\Desktop\GMRTv3_7_20200720topo.grd');

% Pre-processing of time series to calculate wave number (k) and then Wave Energy Flux (P)
g=9.87; %gravity m/s^2
dens=1025; %density kg/m^3
Ne=13;  %number of stations
for i=1:Ne
filename=[filebase,num2str(i),'.nc'];  %Reading each .nc data file

H(i)= ncreadatt(filename,'/','Depth'); % Depth of the continental shelf at each station
H1=abs(H(i));

latw(i) =double( ncreadatt(filename,'/','Latitude')); %Latitude
longw(i) = double(ncreadatt(filename,'/','Longitude'));  %Longitude

T = ncread(filename,'time'); %Time
time = double(T);
time = time(2922:end-1);
times=datetime(time,'ConvertFrom','datenum'); %Convert to date format

hs = ncread(filename,'hs');  % Significant height (Hs)
hs = double(hs)/1000;  % convert mm to m
hs=hs(2922:end-1);
HS(:,i)=hs(~isnan(hs));

tp = ncread(filename,'tp');  % Peak period (Tp)
tp = double(tp)/1000;  % convert ms to s
tp=tp(2922:end-1);
TP(:,i)=tp(~isnan(tp));

dp = ncread(filename,'dp'); % Peak direction, azimuthal degrees
dp = double(dp)/100;
dp=dp(~isnan(dp));
dp = dp(2922:end-1);
dprad = (dp).*pi/180;
DP(:,i)=dp;
DP_rad(:,i)=dprad;

% Calculation of Wave energy flux by Newton-Raphson solution of dispersion eq.  
for ii=1:length(hs)
    Ew(ii)=1/8*dens*g*hs(ii).^2;   % Wave energy
    [k(ii),iter,e_a] = jfpa_dispersionNewtonRaphson(g,H1,tp(ii));  %Dispersion equation solved by Newton-Raphson
    Cg(ii)=sqrt(g*k(ii)*tanh(k(ii)*H1))/k(ii)*1/2*(1+2*k(ii)*H1/sinh(2*k(ii)*H1)); %Wave group celerity
    P(ii)=Ew(ii).*Cg(ii);   %Wave energy flux
    DoC_i(ii)=2.28*hs(ii)-68.5*((hs(ii)^2)/(g*tp(ii)^2));
    DoC_o(ii)=(hs(ii)-0.3*0.2834)*tp(ii)*(g/(5000*0.001))^0.5;
end
PP(:,i)=P;
DoCi(:,i)=DoC_i;
DoCo(:,i)=DoC_o;

%variabilidad oleaje
meanP(i)=sum(P)/length(P);
COV(i)=std(P)/mean(P); %coeficiente de variabilidad
i
end

%% Part 1: Wave Climate variability and statistical trends (Classical description and exploratory analysis)
%% Hs-Tp diagrams
station=1;
nPoints=length(HS(:,station));

X(1,:)=HS(:,station); % Significant height stored in X
X(2,:)=TP(:,station); % Peak period stored in X
Xmean=mean(X,2);  % Averaged time series
B=X-Xmean*ones(1,nPoints);  % Normalization
[U,S,V]=svd(B/sqrt(nPoints),'econ');  % Singular value decomposition
theta2=(0:0.01:1)*2*pi;  % Assigning the elipsoid angle  
Xstd=U*S*[cos(theta2); sin(theta2)]; % Constructing the Standard Deviation of the Principal Component Analysis

figure()
subplot(2,1,1)
scatter(X(1,:),X(2,:),'k.','LineWidth',2)
box on, grid on, hold on
plot(Xmean(1)+Xstd(1,:),Xmean(2)+Xstd(2,:),'r-','LineWidth',1.5)
hold on
plot(Xmean(1)+2*Xstd(1,:),Xmean(2)+2*Xstd(2,:),'r-','LineWidth',1.5)
hold on
plot([Xmean(1) Xmean(1)+U(1,1)*S(1,1)],[Xmean(2) Xmean(2)+U(2,1)*S(1,1)],'c-','LineWidth',2)
hold on
plot([Xmean(1) Xmean(1)+U(1,2)*S(2,2)],[Xmean(2) Xmean(2)+U(2,2)*S(2,2)],'c-','LineWidth',2)
xlabel('Significant height (H_s) [m]','FontName','times')
ylabel('Peak period (T_m) [s]','FontName','times')
subplot(2,1,2)
polarscatter(X(1,:),X(2,:),'k.','LineWidth',2)
%% Scatter 3D

figure()
subplot(2,2,1)
scatter3(HS(TP(:,1)<=10,1),TP(TP(:,1)<=10,1),DP(TP(:,1)<=10,1),'MarkerEdgeColor',[0 .75 .75],'MarkerFaceColor','r')
hold on
scatter3(HS(TP(:,1)>10,1),TP(TP(:,1)>10,1),DP(TP(:,1)>10,1),'MarkerEdgeColor',[0 .75 .75],'MarkerFaceColor','g')
subplot(2,2,2)
scatter3(HS(TP(:,5)<=10,5),TP(TP(:,5)<=10,5),DP(TP(:,5)<=10,5),'MarkerEdgeColor',[0 .75 .75],'MarkerFaceColor','r')
hold on
scatter3(HS(TP(:,5)>10,5),TP(TP(:,5)>10,5),DP(TP(:,5)>10,5),'MarkerEdgeColor',[0 .75 .75],'MarkerFaceColor','g')
subplot(2,2,3)
scatter3(HS(TP(:,8)<=10,8),TP(TP(:,8)<=10,8),DP(TP(:,8)<=10,8),'MarkerEdgeColor',[0 .75 .75],'MarkerFaceColor','r')
hold on
scatter3(HS(TP(:,8)>10,8),TP(TP(:,8)>10,8),DP(TP(:,8)>10,8),'MarkerEdgeColor',[0 .75 .75],'MarkerFaceColor','g')
subplot(2,2,4)
scatter3(HS(TP(:,13)<=10,13),TP(TP(:,13)<=10,13),DP(TP(:,13)<=10,13),'MarkerEdgeColor',[0 .75 .75],'MarkerFaceColor','r')
hold on
scatter3(HS(TP(:,13)>10,13),TP(TP(:,13)>10,13),DP(TP(:,13)>10,13),'MarkerEdgeColor',[0 .75 .75],'MarkerFaceColor','g')
xlabel('H_s (m)')
ylabel('T_p (s)')
zlabel('D_p (°)')
legend('Swell','Local')
%% WaveRose by Ashton & Nienhuis, later by Paniagua-Arroyave
% ref_coast=180;
theta_shore=[130 130 214 214 214 130 180 180 130 120 130 130 130]; %shoreline orientation
% theta_shore=theta_shore;
for i=1:13
[Qw_max(:,i), Qw_net(:,i), Qw_uni(i,:), E(:,i), E_shore(:,i), theta(i,:), phi(i,:)] = jfpa_Qwave(theta_shore(i), DP(:,i), HS(:,i), TP(:,i));
end

st=1; %station
subplot(1,2,2)
DP_stat=DP(:,st);
dp_s=DP_stat-180;
% histogram(E_shore,100,'Normalization','probability','BinLimits',[-180 180],'FaceColor','#00FFFF')
bar(theta(st,:),E_shore(:,st)*10,'FaceColor','k')
xlabel('\phi_{0} (°)')
ylabel('E(%)')
xticks([-90 -45 0 45 90])
subplot(1,2,1)
plot(theta(st,:),Qw_net(:,st),'k','LineWidth',2)
xlabel('\phi_{0} (°)')
ylabel('Q_{s} (kg/s)')
xticks([-90 -45 0 45 90])
% axis([0 360 min()])
%% Latitudinal changes of statistical parameters
[pmund,s]=polyfit(latw,meanP,1);
[y,delta]=polyval(pmund,latw,s);
Md=fitlm(latw,meanP);
R2=Md.Rsquared.Ordinary;  %R cuadrado
p1=fminsearch('fit1',[1 1],[],latw,meanP);
y2=polyval(p1,latw);
sz = linspace(10,100,Ne);

figure()
yyaxis left
scatter(latw,meanP,sz,'b^')
hold on
plot(latw,y2,'k',latw,y2+(10+delta),':r',latw,y2-(10+delta),':b','LineWidth',1);
hold on
plot(latw,meanP,'-b')
set(gca,'YColor','b')
ylabel('Energy flux  (P) [W/m]','FontName','times')
xlabel('Latitude [°]','FontName','times')
yyaxis right
scatter(latw,COV,sz,'r^')
set(gca,'YColor','r')
ylabel('Variability coefficients (COV) [%]','FontName','times')

%% Part 2: Coastal wave energy flux availability regimes in the continental shelf (Classification and description)
%% Principal Component Space (from PCA)
nPoints=length(PP);
X=PP';
Xmean=mean(X,2);
B=X-Xmean*ones(1,nPoints);
[U2,S2,V2]=svd(B/sqrt(nPoints),'econ');
figure(1), hold on
for i=1:size(X,1)
    x(i)=V2(:,1)'*X(i,:)';
    y(i)=V2(:,2)'*X(i,:)';
    z(i)=V2(:,3)'*X(i,:)';
    if i <= 5  % Low latitudes dissipation regime (station 1-5)
        plot3(x(i),y(i),z(i),'r^','LineWidth',3)
    elseif i>5 && i<12  % Mid latitude dissipation regime (station 5-11)
        plot3(x(i),y(i),z(i),'g^','LineWidth',3)
    else   % High latitude dissipation regime (station 12 and 13)
        plot3(x(i),y(i),z(i),'b^','LineWidth',3)
    end
end
grid on, box on
xlabel('PC 1','FontName','times')
ylabel('PC 2','FontName','times')
zlabel('PC 3','FontName','times')
set(gca,'FontSize',15), axis tight

figure(2)
subplot(1,2,1)
semilogy(diag(S2),'k-o','LineWidth',2.5)
set(gca,'FontSize',15), axis tight
box on, grid on
xlabel('PCs','FontName','times')
ylabel('Singular values','FontName','times')
subplot(1,2,2)
plot(cumsum(diag(S2))./sum(diag(S2)),'k-o','LineWidth',2)
set(gca,'FontSize',15), axis tight
box on, grid on
xlabel('PCs','FontName','times')
ylabel('Cumulative % of variance','FontName','times')

%% Naive bayes classification and K-mean clustering
cell={'zone 1','zone 1','zone 1','zone 1','zone 1','zone 2','zone 2','zone 2','zone 2','zone 2','zone 2','zone 3','zone 3'}';
U=cat(2,x',y');
Model = fitcnb(U,cell,'ClassNames',{'zone 1','zone 2','zone 3'});

figure(1)
subplot(1,2,1)
gscatter(U(:,1),U(:,2),cell);
h = gca;
cxlim = h.XLim;
cylim = h.YLim;
hold on
Params = cell2mat(Model.DistributionParameters); 
Mu = Params(2*(1:3)-1,1:2); % Extract the means
Sigma = zeros(2,2,3);
for j = 1:3
    Sigma(:,:,j) = diag(Params(2*j,:)).^2; % Create diagonal covariance matrix
    xlim = Mu(j,1) + 4*[-1 1]*sqrt(Sigma(1,1,j));
    ylim = Mu(j,2) + 4*[-1 1]*sqrt(Sigma(2,2,j));
    f = @(x,y) arrayfun(@(x0,y0) mvnpdf([x0 y0],Mu(j,:),Sigma(:,:,j)),x,y);
    fcontour(f,[xlim ylim]) % Draw contours for the multivariate normal distributions 
end
h.XLim = cxlim;
h.YLim = cylim;
title('Naive Bayes Classifier','FontName','times')
xlabel('PC 1','FontName','times')
ylabel('PC 2','FontName','times')
legend('class 1','class 2','class 3','FontName','times')
hold off
% axis([-8e6 1 -1.5e6 2e6])
% K mean clustering
X(:,1)=U(:,1);
X(:,2)=U(:,2);
opts = statset('Display','final');
[idx,C] = kmeans(X,3,'Distance','cityblock',...
    'Replicates',5,'Options',opts);
set(gca,'FontSize',15), axis tight
subplot(1,2,2)
plot(X(idx==1,1),X(idx==1,2),'b*','MarkerSize',12)
hold on
plot(X(idx==2,1),X(idx==2,2),'r*','MarkerSize',12)
hold on
plot(X(idx==3,1),X(idx==3,2),'g*','MarkerSize',12)
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('cluster 1','cluster 2','cluster 3','centroids',...
       'Location','NW')
title('K-means Clustering','FontName','times')
xlabel('PC 1','FontName','times')
ylabel('PC 2','FontName','times')
hold off
set(gca,'FontSize',15), axis tight
% axis([-8e6 1 -1.5e6 2e6])

%% Part 3: Climate periodicities and spectral correlation and coherence to the Oceanic Niño Index (Description)
%% Time series analysis, DMD, Wavelet and Fourier transforms, and Wavelet Coherence
%% Dynamic Mode Decomposition
X=PP';
X1=X(:,1:end-1);
X2=X(:,2:end-1);
r=4;
dt=(time(2)-time(1))*24*60*60; % dt en segundos
ds=50000;  % ds en metros
% max_cyc=12;
% L=4;
[Phi, S, V]=svd(X,'econ');
[Phi1, omega1, lambda1, b, b2, Xdmd, timedy,f,P] = DMD(X1,X2,r,dt);
PP_rec=Xdmd';

% [tree, omega1, b]=mrDMD(X, dt, r, max_cyc, L);
% [map, low_f] = mrDMD_map(tree);

dominantMode=find(b==max(b));

% figure(1)
% stem(f, P, 'k');
% set(gca, 'XScale', 'log')
% title('DMD Spectrum')
% xlabel('f [Hz]');ylabel('Amplitude')
% subplot(2,1,1)
% plot(Vr)
% subplot(2,1,2)
% plot(PP_rec(:,1:4))

figure('color', 'w');
stem(abs(omega1),abs(b), 'k^', 'filled'); hold on
plot(abs(omega1(dominantMode)),abs(b(dominantMode)),'ro','MarkerSize',10)
set(gca, 'XScale', 'log')
title('DMD Spectrum')
xlabel('f [Hz]');ylabel('Amplitude')
%% blocks fourier
dts=60*60*3/(60*60*24*30*12); % 3 horas en segundos para frecuencia de sampleo
% dts=3/360;
MSLR=detrend(PP(:,13));

elem=2^11;
solap=50;
M=floor(length(MSLR)/elem)/(1-solap/100)-1;
clear indices

for jj=1:M
    index1=(jj-1)*elem*(1-solap/100)+1;
    index2=(jj-1)*elem*(1-solap/100)+elem;
    indices(jj,1)=index1;
    indices(jj,2)=index2;
% % %     jj
end

f_blo=(1:1:elem)'/(elem*dts);
N=length(elem);
S=nan(elem/2,M);
for ii=1:M
    y_b1=MSLR((indices(ii,1):indices(ii,2)));
    for iii=1:length(f_blo)
        Ytest(:,iii)=sqrt(8/3)*dts*sum(y_b1.*exp(-1i*2*pi*f_blo(iii)*(1:1:length(y_b1))*dts)');
    end
    S(:,ii)=2/(N*dts)*Ytest(1:elem/2).*conj(Ytest(1:elem/2));
% % %     ii
end

totb=M+(M-1)*(1/(1-solap/100)-1);
dof=totb*2*2;
conf=0.95;
conf_int=[dof/chi2inv(conf,dof) dof/chi2inv(1-conf,dof)];

figure(1)
plot(f_blo(1:elem/2),mean(S,2),'k.-')
hold on
plot([1080 1080],1.107e7*conf_int,'-r','LineWidth',2)
hold on
plot([360 360],4.7e8*conf_int,'-r','LineWidth',2)
hold on
plot([1.406 1.406],1.17e10*conf_int,'-r','LineWidth',2)
% grid on
xlabel('Frequency (cpy)')
ylabel('Spectral density (W^2/m^2/cycles/year)')
set(gca,'xscale','log','yscale','log')
% grid on
% box on
%% Time series analysis
station=1;
X=PP';  % Esta es la matriz de datos
[U2,S2,V2]=svd(X,'econ');  %PCA modes

% seasonal filter
dt=60*60*3; % 3 horas en segundos para frecuencia de sampleo
fs=1/dt;
fm=1/(30*24*60*60);   % 1 mes en Hz
fa=1/(12*30*24*60*60); % 1 año en Hz

[b1,a]=butter(1,fm/(fs/2),'low');
SV=filtfilt(b1,a,PP(:,station));   % Seasonal variability
SV=SV';
Modes=filtfilt(b1,a,V2);   % Low-pass filter of the V modes

% remove outliers
out=[0 97];
[OV_out,TF] = rmoutliers(PP(:,station),'percentiles',out); % Ordinary variability
for pp=1:length(PP(:,station))
if TF(pp)==1
   time_out(pp)=nan; 
else
    time_out(pp)=time(pp);
end
end
time_out=time_out(~isnan(time_out));
dts=min(diff(time));
t=(time(1):dts:time(end));
OV=interp1(time_out,OV_out,t,'linear');  % Interpolation of Ordinary Variability




%% Fourier Analysis
station=1;
P_p=PP(:,station);
Pp(1,:)=P_p(1:20000);
Pp(2,:)=P_p(20000:40000-1);
Pp(3,:)=P_p(40000:60000-1);
Pp(4,:)=P_p(60000:80000-1);
for i=1:4
N=length(Pp(i,:));
dt=(time(2)-time(1))*24*60*60; %sampling period in seg
P_detrend=normalize(Pp(i,:));
fhat=fft(P_detrend,N);
PSD(:,i)=fhat.*conj(fhat)/N;
freq(:,i)=1/(dt*N)*(0:N);
L=1:floor(N/2);
end
freqt=mean(freq,2);
psd=mean(PSD,2);
%% normalization
MMN=(SV-mean(SV()));
n=243.5;
time_new=arrayfun(@(i) mean(time(i:i+n-1)),1:n:length(time)-n+1)'; % the averaged vector
MMN_new=arrayfun(@(i) mean(MMN(i:i+n-1)),1:n:length(MMN)-n+1)'; % the averaged vector
OV_new=arrayfun(@(i) mean(OV(i:i+n-1)),1:n:length(OV)-n+1)'; % the averaged vector
SV_new=arrayfun(@(i) mean(SV(i:i+n-1)),1:n:length(SV)-n+1)'; % the averaged vector
Mode=arrayfun(@(i) mean(Modes(i:i+n-1,:)),1:n:length(Modes)-n+1)'; % the averaged vector
Mode=cell2mat(Mode);
% MMN_new(:,1)=boxpdf(MMN_new(:,1));
t = (1:1:length(Mode))';
Lim_X = [min(t) max(t)];  % Wavelet plot years
Tick_X = (min(t):48:max(t));
TickLab_X = {'Jan-80';'Jan-84';'Jan-88';'Jan-92';'Jan-96';'Jan-00';'Jan-04';'Jan-08'};

%% Figures
figure(1)
subplot(1,2,1)
plot(freqt(L)*(60*60*24*30*12),psd(L),'k','LineWidth',3)
set(gca,'XScale','log','YScale','log')
% hold on
% semilogx(freqt(L),PSDclean(L),'r','LineWidth',3)
legend('raw','clean')
title('Fourier Spectrum')
xlabel('frequency (cpy)');ylabel('Spectral density (W^2/m^2/cycles/year)')
subplot(1,2,2)
stem(abs(omega1)*(60*60*24*30*12),abs(b), 'k^', 'filled'); hold on
plot(abs(omega1(dominantMode)*(60*60*24*30*12)),abs(b(dominantMode)),'ro','MarkerSize',10);
set(gca, 'XScale', 'log')
title('DMD Spectrum')
xlabel('frequency (cpy)');ylabel('Amplitude (W/m)')

figure(3)
subplot(3,1,1)
plot(time,PP(:,station),'k','LineWidth',2)
hold on
plot(time,OV,'r','LineWidth',2)
hold on
plot(time,SV,'g','LineWidth',2)
legend('Raw','OV','SV')
subplot(3,1,2)
plot(time,MMN,'k','LineWidth',2)
legend('MMN')
subplot(3,1,3)
plot(time,Modes(:,1),'b','LineWidth',2)
hold on
plot(time,Modes(:,2),'y','LineWidth',2)
hold on
plot(time,Modes(:,3),'k','LineWidth',2)
hold on
plot(time,Modes(:,4),'c','LineWidth',2)
legend('PC1','PC2','PC3','PC4')


