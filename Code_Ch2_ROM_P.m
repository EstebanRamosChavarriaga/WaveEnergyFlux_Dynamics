%% WAVE ENERGY FLUX VARIABILITY ALONG CONTINENTAL SHELVES FROM WAVE CLIMATE DATA-DRIVEN ANALYTICS. RAMOS-CHAVARRIAGA ET AL., 2021. (Esteban Ramos-Chavarriaga 2021)
% Part 1: Data-driven discovery of nonlinear terms based on spatial
% derivates of the system. (time series analysis and sparse regression)
% Part 2: Reduced-order dynamics of wave energy flux variability along the Northern Andes Pacific coast (Description & Prediction)
%% Pre-processing: RUN THIS SECTION BEFORE RUNING EITHER PART 1 OR 2
% Data acquisition
clc, clear all
% Datos WW3
filebase='C:\Users\eramo\Desktop\Esteban\EAFIT\Maestría en Ciencias de la Tierra\1. Source to Sink\Wave Climate and Energy Dissipation\Data\Data_wavewatch_Pacifico\pacifico_sur_nort_';
% filebase='C:\Users\eramo\Desktop\Esteban\EAFIT\Maestría en Ciencias de la Tierra\1. Source to Sink\Wave Climate and Energy Dissipation\Proyecto final\Data_wavewatch_Pacifico\New_data\WaveWatch_timeseries_colombia_shelf';%% Pre-processing
g=9.87; %gravity m/s^2
dens=1025; %density kg/m^3  
Ne=13;  %number of stations
for i=1:Ne
filename=[filebase,num2str(i),'.nc'];  %Reading each .nc data file

H = ncreadatt(filename,'/','Depth'); % Depth of the continental shelf at each station
H1=abs(H);

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
DP(:,i)=dprad;
%% Calculation of Wave energy flux by Newton-Raphson solution of dispersion eq.  
for ii=1:length(hs)
    Ew(ii)=1/8*dens*g*hs(ii).^2;   % Wave energy
    [k(ii),iter,e_a] = jfpa_dispersionNewtonRaphson(g,H1,tp(ii));  %Dispersion equation solved by Newton-Raphson
    Cg(ii)=sqrt(g*k(ii)*tanh(k(ii)*H1))/k(ii)*1/2*(1+2*k(ii)*H1/sinh(2*k(ii)*H1)); %Wave group celerity
    P(ii)=Ew(ii).*Cg(ii);   %Wave energy flux
end
PP(:,i)=P;
kappa(:,i)=k;
%variabilidad oleaje
meanP(i)=sum(P)/length(P);
COV(i)=std(P)/mean(P); %coeficiente de variabilidad
% meanPP(i) = meanP(meanP~=0);
% COVP(i) = COV(COV~=0);
i
sz = linspace(10,100,Ne);
end
%% PART 1: DATA-DRIVEN DISCOVERY OF NONLINEAR DYNAMICS
% Always run Filters before any other time series computation.
%% Time series analysis: Obtaining different training datasets to discover and solve the model
%% Filters
T=60*60*3; % 3 horas en segundos para frecuencia de sampleo
fs=1/T; 
fm=1/(30*24*60*60);   % 1 mes en Hz
fa=1/(12*30*24*60*60);   % 1 año en Hz
fd=1/(5*12*30*24*60*60);
[b3,a3]=butter(1,fa/(fs/2),'low');
PP_low=filtfilt(b3,a3,PP);
PP_high=PP;
nn=360;
n=length(PP)/nn;
for pp=1:Ne
PPP=PP_low(:,pp);
PPPP=PP_high(:,pp);
PP_filt(:,pp)=arrayfun(@(i) mean(PPP(i:i+n-1)),1:n:length(PPP)-n+1)'; % the averaged vector
PP_new_ordinary(:,pp)=arrayfun(@(i) mean(PPPP(i:i+n-1)),1:n:length(PPPP)-n+1)'; % the averaged vector
end
time_low=arrayfun(@(i) mean(time(i:i+n-1)),1:n:length(time)-n+1)'; % the averaged vector
%% Random sample
nn=360;
y=randsample(length(PP_filt(:,1)),nn);
for q=1:Ne
P_sel=PP_filt(:,q);
PP_rand(:,q)=P_sel(y);
end
n=length(time)/nn;
time_low=arrayfun(@(i) mean(time(i:i+n-1)),1:n:length(time)-n+1)'; % the averaged vector
%% PSC
% denoise data fft Fourier Analysis
for station=1:Ne
N=length(PP_filt(:,station));
dt=(time(2)-time(1))*24*60*60; %sampling period in seg
P_detrend=normalize(PP_filt(:,station));
fhat=fft(P_detrend,N);
PSD=fhat.*conj(fhat)/N;
freq=1/(dt*N)*(0:N);
L=1:floor(N/2);
indices=PSD>(mean(PSD)+std(PSD))/0.2;
% indices=PSD>100;
PSDclean=PSD.*indices;
fhat=indices.*fhat;
PSC(station,:)=ifft(fhat);
end
nn=360;
n=length(PP)/nn;
time_low=arrayfun(@(i) mean(time(i:i+n-1)),1:n:length(time)-n+1)'; % the averaged vector
%% Uniform sample
y1=zeros(1,length(PP_filt(:,1)));
y1(1:10)=0;
y1(10:300)=1;
y1(300:length(PP_filt(:,1)))=0;
for qq=1:Ne
for qqq=1:length(y1)
if y1(qqq)==0
PP_uni(qqq,qq)=nan;
elseif y1(qqq)==1
PP_uni(qqq,qq)=PP(qqq,qq);
end
end
end
for i=1:Ne
Pp=PP_uni(:,i);
PP_unif(:,i)=Pp(~isnan(Pp));
end
%% Model discovery PDE-FIND algorithmv(spatial derivates and sparse regression)
%% Numerical differentiation schemes
% P_new=cat(2,PP(:,11),PP(:,25),PP(:,41),PP(:,5),PP(:,68),PP(:,81),PP(:,90),PP(:,97)); %Initial system one dimension 30° arriving waves
% P_new=PP';
% P_new=PP_filt';
% P_new=PP_new_ordinary';
% P_new=PP_rand;
% P_new=PSC;
P_new=PP_unif;

[m,n]=size(P_new);
dt=(time(2)-time(1))*24*60*60; % dt en segundos

P_t=zeros(m,n-2);
for jj=1:m
   for j=2:n-1
      P_t(jj,j-1)=(P_new(jj,j+1)-P_new(jj,j-1))/(2*dt); % temporal derivates
   end
end

dx=50000;

D=zeros(m,m); D2=zeros(m,m);
for j=1:m-1
   D(j,j+1)=1;
   D(j+1,j)=-1;
   
   D2(j,j+1)=1;
   D2(j+1,j)=1;
   D2(j,j)=-2;
end
D(m,1)=1;
D(1,m)=-1;
D=(1/(2*dx))*D;

D2(m,m)=-2;
D2(m,1)=1;
D2(1,m)=1;
D2=D2/(dx^2);

for jj=2:n-1
   P_x(jj-1,:)=((D*P_new(:,jj)));
   P_xx(jj-1,:)=((D2*P_new(:,jj)));
   P2_x(jj-1,:)=((D*(P_new(:,jj).^2)));
   jj
end

P_0=reshape(P_new(:,2:end-1),(n-2)*m,1); % first column of the A matrix
Pt=reshape(P_t,(n-2)*m,1); 
Px=reshape((P_x'),(n-2)*m,1);
Pxx=reshape((P_xx'),(n-2)*m,1); 
P2x=reshape((P2_x'),(n-2)*m,1);

%% Solution of Ax=b to obtain library parameters
A=[P_0 P_0.^2 P_0.^3 Px Px.*Px Pxx Pxx.*Pxx P2x P_0.*Px P_0.*Pxx cos(P_0) cos(Px) cos(Pxx) cos(P2x) cos(Px.*Px) cos(Pxx.*Pxx)];
% A=[P_0 P_0.^2 P_0.^3 Px Pxx Pxx P2x P_0.*Px P_0.*Pxx Px.*Px];
% A=[Px Pxx P2x Px.*P_0 Px.*Px Px.*Pxx];
% A=[Px Px.*Px P2x cos(Px) sin(Px) sin(Px.*Px) cos(Px.*Px)];
% A=[P_0 P_0.^2 Px P2x P_0.*Px];
x1=A\Pt;
% x1=lasso(A,Pt,'Alpha',0.01,'Lambda',0.0001);
figure()
barh(x1,'FaceColor','#00FFFF','EdgeColor','k','LineWidth',1.5)
% discovered equation =0.6*(phi')*phi_x*x+0.8*(phi')*(phi_x*x).^2+0.01*(phi')*cos((phi_x*x).^2)
%% PART 2: REDUCED-ORDER MODEL (ROM) OF THE SYSTEM (RUN THIS PART ONLY WITH ONE (1) TIME SERIES)
%% Principal component analysis (PCA) and dynamic mode decomposition (DMD) to obtain the optimal basis set
%% Singular Value Decomposition
% [U,S,V]=svd(PP','econ');
% figure(2)
% plot(diag(S)/sum(diag(S)),'ko','Linewidth',2)
% 
% figure(3)
% plot(U(:,1:3),'Linewidth',2)
% figure(4)
% plot(V(:,1:3),'Linewidth',2)
%% Dynamic Mode Decomposition
%original
% PP_train=PP(1:end-28800,:);  % data from 1980 to 2000
% PP_test=PP(end-28799:end,:);  % data from 2000 to 2010
% t_train=time(1:end-28800,:);
% t=time(end-28799:end,:);

%spectral components
% PP_train=PSC(:,1:end-abs(nn/4))';  % data from 1980 to 2000
% PP_test=PSC(:,1:end-abs(nn/4))';  % data from 1980 to 2000
% t_train=time(:,1:end-abs(nn/4))';  % data from 1980 to 2000
% t=time_low(end-abs(nn/4):end,:);

%seasonal
PP_train=PP_filt(1:end-abs(nn/4),:);  % data from 1980 to 2000
PP_test=PP_filt(end-abs(nn/4):end,:);  % data from 2000 to 2010
t_train=time_low(1:end-abs(nn/4),:);
t=time_low(end-abs(nn/4):end,:);
% PP_train=PP_new_ordinary(1:end-abs(nn/4),:);  % data from 1980 to 2000
% PP_test=PP_new_ordinary(end-abs(nn/4):end,:);  % data from 2000 to 2010
% t_train=time_low(1:end-abs(nn/4),:);
% t=time_low(end-abs(nn/4):end,:);

%random
% PP_train=PP_rand(1:end-abs(nn/4),:);  % data from 1980 to 2000
% PP_test=PP_rand(end-abs(nn/4):end,:);  % data from 2000 to 2010
% t_train=time_low(1:end-abs(nn/4),:);
% t=time_low(end-abs(nn/4):end,:);

%uniform
% PP_train=PP_unif(1:end-abs(nn/4),:);  % data from 1980 to 2000
% PP_test=PP_unif(end-abs(nn/4):end,:);  % data from 2000 to 2010
% t_train=time_low(1:end-abs(nn/4),:);
% t=time_low(end-abs(nn/4):end,:);

% %seasonal future
% PP_train=PP_rand;  % data from 1980 to 2000
% t=(time_low(end):time_low(2)-time_low(1):time_low(end)+100000);
% times=datetime(t,'ConvertFrom','datenum'); %Convert to date format

X = PP_train';
X1 = X(:,1:end-1);
X2 = X(:,2:end-1);
r=5;
dt=(time(2)-time(1))*24*60*60; % dt en segundos
[U_d,S_d,V_d]=svd(X1,'econ');
Ur=U_d(:,1:r);
Sr=S_d(1:r,1:r);
Vr=V_d(:,1:r);
Vr=Vr(1:end-1,:);
Atilde=Ur'*X2*Vr/Sr;
[W,D]=eig(Atilde);
Phi=X2*Vr/Sr*W; %dmd modes
lamda=diag(D);
omega=log(lamda)/dt;
%% plots
figure(1)
subplot(2,1,1); plot(real(Vr),'LineWidth',3);
legend('mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5')
set(gca,'FontSize',15), axis tight
box on, grid on
xlabel('Coastal stations','FontName','times')
ylabel('Normalized Amplitude','FontName','times')
subplot(2,1,2); plot(U_d(:,1:r),'LineWidth',3);
legend('mode 1', 'mode 2', 'mode 3', 'mode 4', 'mode 5')

figure(2)
subplot(1,2,1)
plot(cumsum(diag(S_d))./sum(diag(S_d)),'k-o','LineWidth',2)
set(gca,'FontSize',15), axis tight
box on, grid on
xlabel('Principal components (PCs)','FontName','times')
ylabel('Cumulative percent of variance','FontName','times')
subplot(1,2,2)
semilogy(diag(S_d),'k-o','LineWidth',2.5)
set(gca,'FontSize',15), axis tight
box on, grid on
xlabel('Principal components (PCs)','FontName','times')
ylabel('Singular values','FontName','times')

figure(3)
hold on
for i=1:size(X1,1)
    x(i)=V_d(:,1)'*X1(i,:)';
    y(i)=V_d(:,2)'*X1(i,:)';
    z(i)=V_d(:,3)'*X1(i,:)';
    if i<=5
        plot3(x(i),y(i),z(i),'bo','LineWidth',3)
    elseif i>5 && i<12
        plot3(x(i),y(i),z(i),'ro','LineWidth',3)
    elseif i>11 && i<=13
        plot3(x(i),y(i),z(i),'co','LineWidth',3)
    elseif i==21
        plot3(x(i),y(i),z(i),'yo','LineWidth',3)
    else
        plot3(x(i),y(i),z(i),'ko','LineWidth',3)
    end
end
grid on, box on
xlabel('PC 1','FontName','times')
ylabel('PC 2','FontName','times')
zlabel('PC 3','FontName','times')
set(gca,'FontSize',15), axis tight
%% Solution of reduced dynamical system
%% Differentiation shcemes applied to phi_r
% P0=PP(end,:);
% P0=mean(PP);
% P0=([7238.486 11429.196 7795.363 6093.334 1815.189 2219.937 2097.254 5872.103 5784.656 5127.642 6555.744 9680.730 10149.171]);   % Inital contidions 1980
% P0=([0.6315    1.0946    0.9145    1.1056    0.9746    0.9969    1.0455    1.1615    1.2195    1.0620    1.1690    1.8737    2.6998]).*1e+04;   % Inital contidions 2000
P0=PP_train(end,:);
% kappa1=kappa'; 
% kappa2=kappa1(:,1:3);
r=5;
%r=4;  % rank of projection
% phi=PP(1:r,:)';
% phi=real(Phi(:,1:r));  % Phi_r DMD modes
phi=U_d(:,1:r);  % Phi_r POD modes
[m,n]=size(phi);
for j=1:r
%   phi_x(:,j)=real(ifft((i*kappa2(:,j)).*fft(phi(:,j)))); % second derivatives
  a0(j)= P0*conj(phi(:,j));  % projection of initial cond   itions
end

%Spatial derivatives
dx=50000;
D=zeros(m,m); D2=zeros(m,m);
for j=1:m-1
   D(j,j+1)=1;
   D(j+1,j)=-1;
   
   D2(j,j+1)=1;
   D2(j+1,j)=1;
   D2(j,j)=-2;
end
D(m,1)=1;
D(1,m)=-1;
D=(1/(2*dx))*D;

D2(m,m)=-2;
D2(m,1)=1;
D2(1,m)=1;
D2=D2/(dx^2);

for jj=2:n-1
   Phi_x(jj-1,:)=(D*phi(:,jj));
   Phi_xx(jj-1,:)=(D2*phi(:,jj));
   Phi2_x(jj-1,:)=(D*(phi(:,jj).^2));
   jj
end
phi_x=Phi_x';
phi_xx=Phi_xx';
phi_2x=Phi2_x';
%% Numerical Solution
% [t,asol1]=ode45('rhs_model_full',t,a0,[],phi(:,1),phi_x(:,1),phi_xx(:,1),phi_2x(:,1));

[t,asol]=ode45('rhs_model_dummy',t,a0,[],phi,phi_x,phi_xx);
% [t,asol2]=ode45('rhs_model_dummy',t,a0,[],phi(:,2),phi_x(:,2),phi_xx(:,2),phi_2x(:,2));
% [t,asol3]=ode45('rhs_model_dummy',t,a0,[],phi(:,3),phi_x(:,3),phi_xx(:,3),phi_2x(:,3));

% [tt,asol]=ode45(@(tt,x)ch_rom_a_rhs_burger(tt,x,phi,phi_x,phi_2x),t,a0);
% [t,asol]=ode45('ch_rom_a_rhs_less',t,a0,[],phi,phi_x);  %discovered 2
% [t,asol]=ode45('ch_rom_a_rhs_more',t,a0,[],phi,phi_x);  %discovered 6
% [t,asol]=ode45('ch_rom_a_rhs',t,a0,[],phi,phi_2x);       %discovered 4
% [t,asol]=ode45('ch_rom_a_rhs_low',t,a0,[],phi,phi_x);       %discovered low
% us=phi(:,1)*asol(:,1)'+phi(:,2)*asol(:,2)'+phi(:,3)*asol(:,3)'+phi(:,4)*asol(:,4)';
% us=phi*asol1'+phi*asol2'+phi*asol3';
%
us=phi*asol';
us=us';
% us=phi(:,1)*asol(:,1)';

% figure(1)
% subplot(2,1,1)
% surfl(PP_test), shading interp, colormap(summer);
% xlabel('Nearshore stations','FontName','times')
% ylabel('Time','FontName','times')
% zlabel('Wave Energy Flux [W/m]','FontName','times')
% set(gca,'FontSize',15)
% % axis([0 15 0 2e4 0 5e4])
% subplot(2,1,2)
% surfl(us), shading interp, colormap(summer);
% xlabel('Nearshore stations','FontName','times')
% ylabel('Time','FontName','times')
% zlabel('Wave Energy Flux [W/m]','FontName','times')
% set(gca,'FontSize',15)
% axis([0 15 0 2e4 0 5e4])
% set(gca,'Xdir','reverse')
% waterfall(latw,t,(us.')); colormap([0 0 0]);
% set(gca,'Xlim',[-15 15],'Ylim',[0 2*pi],'Zlim',[0 4],'Xtick',[-15 0 15],'Ytick',[0 pi 2*pi],'Yticklabel',{'0','\pi','2\pi'},'Ztick',[0 2 4],'Fontsize',[15])
%
surf(PP_test,'FaceColor','g');
xlabel('Nearshore stations')
ylabel('Time (months)')
zlabel('Wave Energy Flux (W/m)')
set(gca,'FontSize',15)
% datetick('y','yy','keepticks')
% axis([0 15 0 2e4 0 5e4])
hold on
surf(us,'FaceColor','b');
xlabel('Nearshore stations')
ylabel('Time (months)')
zlabel('Wave Energy Flux (W/m)')
set(gca,'FontSize',15)
legend('Real','Simulated')
%% Error and model evaluation
% A1=normalize(us(end,:));
% Prod1=normalize(PP_test(end,:));
A1=us(end,:);
Prod1=PP_test(end,:);
[pmund,s]=polyfit(A1,Prod1,1);
[y,delta]=polyval(pmund,A1,s);
[a, yest, sst, rr] = xyplot1(A1,Prod1);
Md=fitlm(A1,Prod1);
R22=Md.Rsquared.Ordinary  %R cuadrado
% RE=sum(Md.Residuals.Standardized)*100
% Error=abs(sum(abs(normalize(Prod1)-normalize(A1))/normalize(Prod1))*100)
% Error=sum(Md.Residuals.Pearson)*100
p1=fminsearch('fit1',[1 1],[],A1,Prod1);
p2=fminsearch('fit2',[1 1],[],A1,Prod1);
p3=fminsearch('fit3',[1 1],[],A1,Prod1);
y3=polyval(p3,A1);
% y3=polyval(p2,A1);
% y3=polyval(p3,A1);
% cov=100-std(y3)/mean(y3)*100
error=abs(sqrt((1/(length(A1)))*sqrt((normalize(A1.^2)-normalize(Prod1.^2)))))*100
etotal=sum(error)/length(error)
error1=Error1(A1, Prod1,'Target and Predicted');


figure(3)
subplot(2,1,1)
scatter(A1,Prod1,'k^','LineWidth',3);
grid on
hold on
plot(A1,y3,'y',A1,y3+(10+delta),':r',A1,y3-(10+delta),':r','LineWidth',3);
xlabel('Simulated (Model)','FontName','times')
ylabel('Real (Data)','FontName','times')
legend('Data','Regression','95% confidence')
set(gcf,'color','w');
set(gca,'FontSize',15), axis tight
subplot(2,1,2)
plot(error,'.-r','LineWidth',3)
tx=['Error =' num2str(etotal)];
text(4,14,tx,'Color','red','FontSize',15)
ylabel('Error (%)','FontName','times')
xlabel('Coastal stations','FontName','times')
set(gca,'FontSize',15), axis tight
%% Errors
clc, clear all, close all
numb_rand=[100 500 500 1000 1000 5000 10000 20000];
numb=[100 500 1000 5000];

rand_errors=[0.0715 8.7148  1.5270 2.8510 3.8951 1.0740 1.2359 3.9826];
rand_rsquared=[0.7692  0.9394 0.5712 0.8612 0.6785 0.9241 0.8391 0.5723];
rand_MSE=[98418285.71 128994506.46 48535352.87 620531539.12 470491013.18 279736803.46 153063176.09 90645623.20];

filt_errors=[1.1286 1.1181 0.9811 0.9191];
filt_rsquared=[0.9314 0.9554 0.9585 0.9597];
filt_MSE=[158232133.24 114407494.07 113294953.12 112626537.89]./numb;

ord_errors=[0.6340 18.7981 5.7243 35.5316];
ord_rsquared=[0.9530 0.5846 0.6745 0.1004];
ord_MSE=[131448980.04 44465398.38 80830411.53 17954942.37]./numb;
filt_cov=[];

figure(1)
scatter(numb_rand,rand_rsquared,'r*','LineWidth',3)
hold on
scatter(numb,filt_rsquared,'b*','LineWidth',3)
hold on
scatter(numb,ord_rsquared,'g*','LineWidth',3)
ylabel('R^{2}')
xlabel('Number of training samples')
set(gca,'FontSize',15), axis tight
set(gca, 'XScale', 'log')
legend('Random','Seasonal','Ordinary')
set(gcf,'color','w');

% figure(1)
% subplot(3,1,1)
% scatter(numb_rand,rand_rsquared,'r*','LineWidth',3)
% hold on
% scatter(numb,filt_rsquared,'b*','LineWidth',3)
% hold on
% scatter(numb,ord_rsquared,'g*','LineWidth',3)
% ylabel('R_{2}')
% set(gca,'FontSize',15), axis tight
% set(gca, 'XScale', 'log')
% subplot(3,1,2)
% scatter(numb_rand,rand_MSE,'r*','LineWidth',3)
% hold on
% scatter(numb,filt_MSE,'b*','LineWidth',3)
% hold on
% scatter(numb,ord_MSE,'g*','LineWidth',3)
% ylabel('MSE')
% set(gca,'FontSize',15), axis tight
% set(gca, 'XScale', 'log')
% subplot(3,1,3)
% scatter(numb_rand,rand_errors,'r*','LineWidth',3)
% hold on
% scatter(numb,filt_errors,'b*','LineWidth',3)
% hold on
% scatter(numb,ord_errors,'g*','LineWidth',3)
% xlabel('Training data samples')
% ylabel('Error (%)')
% % zlabel('R^{2}')
% legend('Random','Seasonal','Ordinary')
% set(gcf,'color','w');
% set(gca,'FontSize',15), axis tight
% set(gca, 'XScale', 'log')


