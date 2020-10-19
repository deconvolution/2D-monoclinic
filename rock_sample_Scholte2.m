%% input
close all;
clear all;

% load image
A=imread('NS.png');
A=im2double(A)';

% add 10 layers around the original image for air
A=[repmat(A(1,:),[10,1]);
    A;
    repmat(A(end,:),[10,1])];
A=[repmat(A(:,1),[1,10]),A,repmat(A(:,end),[1,10])];

% find solid locations
solid=find(A==1);

% find pore locations
fluid=find(A==0);

% dimensions
dt=10^-9;
dx=10^-5;
dz=10^-5;
nt=4*10^5;
nx=size(A,1);
nz=size(A,2);

% receiver locations [m]
r1=(2:10:(nx-1))*dx;
r3=(ones(size(r1)))*(nz-11)*dz;

% PML layers
lp=10;
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
mu=2.94*10^3*(1490)^2;
lambda=2.94*10^3*(2860)^2-2*mu;

% assign solid with its stiffness
C11(solid)=lambda;
C13(solid)=lambda;
C15(solid)=0;
C33(solid)=lambda+2*mu;
C35(solid)=0;
C55(solid)=mu;
rho(solid)=2.94*10^3;

% assign air with its stiffness
C11(fluid)=1145*340^2;
C13(fluid)=1145*340^2;
C15(fluid)=0;
C33(fluid)=1145*340^2;
C35(fluid)=0;
C55(fluid)=0;
rho(fluid)=1145;

% find surrounding layers
[a1,b1]=meshgrid(1:10,1:nz);
[a2,b2]=meshgrid(nx-10+1:nx,1:nz);
[a3,b3]=meshgrid(1:nx,1:10);
[a4,b4]=meshgrid(1:nx,nz-10+1:nz);
IND_air_layer=sub2ind([nx,nz],[reshape(a1,[1,10*nz]),reshape(a2,[1,10*nz]),reshape(a3,[1,10*nx]),reshape(a4,[1,10*nx])],...
    [reshape(b1,[1,10*nz]),reshape(b2,[1,10*nz]),reshape(b3,[1,10*nx]),reshape(b4,[1,10*nx])]);

% assign surrounding layers with air
C11(IND_air_layer)=1145*340^2;
C13(IND_air_layer)=1145*340^2;
C15(IND_air_layer)=0;
C33(IND_air_layer)=1145*340^2;
C35(IND_air_layer)=0;
C55(IND_air_layer)=0;
rho(IND_air_layer)=1145;

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
s1=fix(nx/2)*dx;
s3=11*dz;

% source frequency [Hz]
freq=5*10^6;

% source signal
singles=rickerWave(freq,dt,nt,M);

% give source signal to x direction
src1=zeros(nt,1);
src1=1*[singles];

% give source signal to z direction
src3=src1;
src3=1*[singles];

% source type. 'D' for directional source. 'P' for P-source.
source_type='P';

% point interval in time steps
plot_interval=500;

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
%% plot recordings
% Choose numer of reiceivers to plot. Nr is the array of receiver number.
Nr=[3,4,5];

figure('name','v3 [m/s]');
for i=1:length(Nr)
    subplot(length(Nr),1,i)
    plot((1:nt)*dt,R3(:,Nr(i)));
    xlabel('t [s]');
    ylabel('v3 [m/s]');
    title(num2str(Nr(i)));
end