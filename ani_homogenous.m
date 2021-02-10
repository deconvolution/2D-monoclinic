%% input
close all;
clear all;

% load image
A=ones(200,200);

n=numel(A);

% add 10 layers around the original image for air
A=[repmat(A(1,:),[12,1]);
    A;
    repmat(A(end,:),[12,1])];
A=[repmat(A(:,1),[1,12]),A,repmat(A(:,end),[1,12])];

% find solid locations
solid=find(A==1);

% find pore locations
fluid=find(A==0);

% dimensions
dt=10^-3;
dx=10;
dz=10;
nt=3000;
nx=size(A,1);
nz=size(A,2);

% PML layers
lp=20;

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
C11(:)=lambda+2*mu;
C13(:)=lambda;
C15(:)=0*lambda;
C33(:)=lambda+2*mu;
C35(:)=0*lambda;
C55(:)=mu;
rho(:)=10^3;
%{
% assign air with its stiffness
C11(fluid)=lambda;
C13(fluid)=lambda;
C15(fluid)=0;
C33(fluid)=lambda;
C35(fluid)=0;
C55(fluid)=0;
rho(fluid)=900;

rho(:,100:end)=2*rho(:,100:end);
%}
% find surrounding layers
[a1,b1]=meshgrid(1:lp+2,1:nz);
[a2,b2]=meshgrid(nx-(lp+2)+1:nx,1:nz);
[a3,b3]=meshgrid(1:nx,1:lp+2);
[a4,b4]=meshgrid(1:nx,nz-(lp+2)+1:nz);
IND_air_layer=sub2ind([nx,nz],[reshape(a1,[1,(lp+2)*nz]),reshape(a2,[1,(lp+2)*nz]),reshape(a3,[1,(lp+2)*nx]),reshape(a4,[1,(lp+2)*nx])],...
    [reshape(b1,[1,(lp+2)*nz]),reshape(b2,[1,(lp+2)*nz]),reshape(b3,[1,(lp+2)*nx]),reshape(b4,[1,(lp+2)*nx])]);
clear a1 a2 a3 a4 b1 b2 b3 b4
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
s1=80;%fix(nx/2);
s3=90;%150;

% source frequency [Hz]
freq=1;

% source signal
singles=rickerWave(freq,dt,nt,M);

% give source signal to x direction
src1=zeros(nt,1);
src1=-1*[singles];

% give source signal to z direction
src3=src1;
src3=1*[singles];

% receiver locations [m]
r1=[25:5:170,25:5:170,ones(size(25:5:170))*25,ones(size(25:5:170))*170];
r3=[ones(size(25:5:170))*25,ones(size(25:5:170))*170,25:5:170,25:5:170];
%r1=[25:5:170];
%r3=[ones(size(25:5:170))*25];

% source type. 'D' for directional source. 'P' for P-source.
source_type='D';

% point interval in time steps
plot_interval=200;

% whether to save figure. 1 for save. 0 for not.
save_figure=1;

% plot source
plot_source=1;

% figure path
p2=mfilename('fullpath');
path=[p2 '/'];

% save wavefield
save_wavefield=1;
%% boundary condition
b1=[];
b3=[];
B1=0;
B3=0;
%% pass parameters to solver
[v1,v3,R1,R3]=monoclinic_2D_xz(dt,dx,dz,nt,nx,nz,...
    r1,r3,...
    s1,s3,src1,src3,source_type,...
    lp,nPML,R,...
    rho,C11,C13,C15,C33,C35,C55,...
    Eta11,Eta13,Eta15,Eta33,Eta35,Eta55,...
    b1,b3, ...
    B1,B3, ...
    plot_interval,plot_source,...
    save_figure,path,...
    save_wavefield);
%% write to gif
sources=path;
delaytime=.2;
filename='animation';
gifmaker(filename,delaytime,sources);
%% plot recordings
%{
% Choose numer of reiceivers to plot. Nr is the array of receiver number.
Nr=[6,10,14,16,18];

figure('name','v3 [m/s]');
for i=1:length(Nr)
    subplot(length(Nr),1,i)
    plot((1:nt)*dt,R3(:,Nr(i)));
    xlabel('t [s]');
    ylabel('v3 [m/s]');
    title(num2str(Nr(i)));
end
%}
%% write recordings
DATA=R1;
save([path 'R1.mat'],'DATA');
save(['./source_reconstruction_ani_homogenous/R1_simu.mat'],'DATA');
DATA=R3;
save(['./source_reconstruction_ani_homogenous/R3_simu.mat'],'DATA');
%%
figure;
subplot(2,2,1)
plot(R1);
subplot(2,2,2)
plot(R3);
subplot(2,2,3)
plot(singles)
subplot(2,2,4)
plot(diff(src3)/dt)
%
max(diff(src1)/dt)/max(diff(src3)/dt)
max(R1)/max(R3)
figure
plot(R3./[diff(src3)/dt;0])
