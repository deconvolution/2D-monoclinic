%% input
close all;
clear all;

% load image
A=ones(150,150);

A(:,1:30)=0;

% find solid locations
solid=find(A==1);

% find pore locations
fluid=find(A==0);

% dimensions
dt=10^-3;
dx=10;
dz=10;
nt=6000;
nx=size(A,1);
nz=size(A,2);

% PML layers
lp=15;

% PML coefficient, usually 2
nPML=2;

% Theoretical coefficient, more PML layers, less R
% Empirical values
% lp=[10,20,30,40]
% R=[.1,.01,.001,.0001]
R=.001;

% generate empty density
rho=zeros(nx,nz);

% generate empty elastic stiffness
C11=zeros(nx,nz);
C13=C11;
C15=C11;
C33=C11;
C35=C11;
C55=C11;

% generate empty viscous modulus
Eta11=C11;
Eta13=C11;
Eta15=C11;
Eta33=C11;
Eta35=C11;
Eta55=C11;

% Lame constants for solid
mu=10^9;
lambda=10^9;

% assign solid with its stiffness
C11(solid)=lambda+2*mu;
C13(solid)=lambda;
C15(solid)=0*lambda;
C33(solid)=lambda+2*mu;
C35(solid)=0*lambda;
C55(solid)=mu;

% assign air with its stiffness
lambdaf=1000*340^2;
C11(fluid)=lambdaf;
C13(fluid)=lambdaf;
C15(fluid)=0;
C33(fluid)=lambdaf;
C35(fluid)=0;
C55(fluid)=0;
rho=1000+sqrt(0)*randn(nx,nz);

% find surrounding layers
[a1,b1]=meshgrid(1:lp+2,1:nz);
[a2,b2]=meshgrid(nx-(lp+2)+1:nx,1:nz);
[a3,b3]=meshgrid(1:nx,1:lp+2);
[a4,b4]=meshgrid(1:nx,nz-(lp+2)+1:nz);
IND_air_layer=sub2ind([nx,nz],[reshape(a1,[1,(lp+2)*nz]),reshape(a2,[1,(lp+2)*nz]),reshape(a3,[1,(lp+2)*nx]),reshape(a4,[1,(lp+2)*nx])],...
    [reshape(b1,[1,(lp+2)*nz]),reshape(b2,[1,(lp+2)*nz]),reshape(b3,[1,(lp+2)*nx]),reshape(b4,[1,(lp+2)*nx])]);

% assign surrounding layers with air
%{
C11(IND_air_layer)=1145*340^2;
C13(IND_air_layer)=1145*340^2;
C15(IND_air_layer)=0;
C33(IND_air_layer)=1145*340^2;
C35(IND_air_layer)=0;
C55(IND_air_layer)=0;
rho(IND_air_layer)=1145;
%}

% ratio for viscous term over stiffness term
scale=0;
Eta11=C11*scale;
Eta13=C13*scale;
Eta15=C15*scale;
Eta33=C33*scale;
Eta35=C35*scale;
Eta55=C55*scale;

% source
% magnitude
M=2.7;

% source locations
rng(1)
random_loc1=randi([lp+1,nx-lp],1,5);
random_loc3=randi([31,nz-lp],1,5);
s1=random_loc1;
s3=random_loc3;

% source frequency [Hz]
freq=10;

% source signal
singles=rickerWave(freq,dt,nt,M);

% give source signal to x direction
src1=repmat(singles,[1,length(s1)]);

% give source signal to z direction
src3=1*repmat(singles,[1,length(s1)]);

% receiver locations [m]
r1=(2:2:(nx-1));
r3=(ones(size(r1)))*31;

% source type. 'D' for directional source. 'P' for P-source.
source_type='P';

% point interval in time steps
plot_interval=100;

% whether to save figure. 1 for save. 0 for not.
save_figure=1;

% figure path
p2=mfilename('fullpath');
path=[p2 '/'];
%% pass parameters to solver
[v1,v3,R1,R3]=monoclinic_2D_xz(dt,dx,dz,nt,nx,nz,...
    r1,r3,...
    s1,s3,src1,src3,source_type,...
    lp,nPML,R,...
    rho,C11,C13,C15,C33,C35,C55,...
    Eta11,Eta13,Eta15,Eta33,Eta35,Eta55,...
    plot_interval,...
    save_figure,path);
%% write to gif
sources=path;
delaytime=.2;
filename='animation';
gifmaker(filename,delaytime,sources);
%% plot recordings
% Choose numer of reiceivers to plot. Nr is the array of receiver number.
Nr=[1:3:15];

figure('name','v3 [m/s]');
for i=1:length(Nr)
    subplot(length(Nr),1,i)
    plot((1:nt)*dt,R3(:,Nr(i)));
    xlabel('t [s]');
    ylabel('v3 [m/s]');
    title(num2str(Nr(i)));
end
%%
ref=70;
int_R1=zeros(size(R1,1)*2-1,size(R1,2));
int_R3=zeros(size(R3,1)*2-1,size(R3,2));
lag1=zeros(1,size(R1,2));
lag3=lag1;

for i=1:size(int_R1,2)
    int_R1(:,i)=xcorr(R1(:,i)',R1(:,ref)');
    int_R3(:,i)=xcorr(R3(:,i),R3(:,ref));
end

figure;
subplot(2,2,1)
imagesc([1,Nr],[dt,dt*nt],R1);
colorbar;
subplot(2,2,2)
imagesc([1,Nr],[-dt*nt,dt*nt],int_R1);
colorbar;
subplot(2,2,3)
imagesc([1,Nr],[dt,dt*nt],R3);
colorbar;
subplot(2,2,4)
imagesc([1,Nr],[-dt*nt,dt*nt],int_R3);
colorbar;
%%
a=zeros(8,8);
a(5:8,1:4)=flip(eye(4),2);
a(5:8,5:8)=eye(4);
b=zeros(15,8);
for i=1:8
    b(:,i)=xcorr(a(:,i),a(:,4));
end
figure;
imagesc(b)