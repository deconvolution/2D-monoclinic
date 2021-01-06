%% input
close all;
clear all;

% load image
A=ones(500,300);

A(:,1:16)=0;

% find solid locations
solid=find(A==1);

% find pore locations
fluid=find(A==0);

% dimensions
dt=10^-3;
dx=5;
dz=5;
nt=10^4;
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
mu=2*10^9;
lambda=10^9;

% assign solid with its stiffness
C11(solid)=lambda+2*mu;
C13(solid)=lambda;
C15(solid)=0*lambda;
C33(solid)=lambda+2*mu;
C35(solid)=0*lambda;
C55(solid)=mu;

% assign air with its stiffness
lambdaf=0;
C11(fluid)=lambdaf;
C13(fluid)=lambdaf;
C15(fluid)=0;
C33(fluid)=lambdaf;
C35(fluid)=0;
C55(fluid)=0;
rng(1)
rho=1000+sqrt(10000)*randn(nx,nz);

rho(:,180:end)=rho(:,180:end)*2/3;
rho(1:200,1:100)=rho(1:200,1:100)*4/3;

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
ns=60000;

% magnitude
M=2*ones(1,ns);

% trigger time
t_trigger=randi([1,nt-2000],[1,ns]);

% sign
si=randi([-1,1],[1,ns]);

% source locations
rng(1)
%{
tt=[20:40:480,20:40:480,ones(size(20:40:280))*50,ones(size(20:40:280))*480];
tt2=[ones(size(20:40:480))*50,ones(size(20:40:480))*280,20:40:280,20:40:280];
%}
s1=randi([20,500-20],[1,ns]);
s3=randi([20,280],[1,ns]);

% source frequency [Hz]
freq=5+sqrt(5)*randn([1,ns]);
freq(freq<1)=1;

% initialize source signal to x direction
src1=zeros(nt,ns);

% initialize source signal to z direction
src3=zeros(nt,ns);

% source signal
for is=1:ns
    src1(t_trigger(is):t_trigger(is)+1999,is)=si(is)*rickerWave(freq(is),dt,2000,M(is));
    src3(t_trigger(is):t_trigger(is)+1999,is)=si(is)*rickerWave(freq(is),dt,2000,M(is));
end

% receiver locations [m]
r1=20:20:480;
r3=(ones(size(r1)))*20;

% source type. 'D' for directional source. 'P' for P-source.
source_type='P';

% point interval in time steps
plot_interval=200;

% whether to save figure. 1 for save. 0 for not.
save_figure=1;

% whether plot source
plot_source=0;

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
    plot_interval,plot_source,...
    save_figure,path);
%% write to gif
sources=path;
delaytime=.2;
filename='animation';
gifmaker(filename,delaytime,sources);
%% plot recordings
% Choose numer of reiceivers to plot. Nr is the array of receiver number.
Nr=[1:10:size(R3,2)];

figure('name','v3 [m/s]');
for i=1:length(Nr)
    subplot(length(Nr),1,i)
    plot((1:nt)*dt,R3(:,Nr(i)));
    xlabel('t [s]');
    ylabel('v3 [m/s]');
    title(num2str(Nr(i)));
end
%%
for ref=1:size(R3,2)
    int_R1=(zeros(size(R1,1)*2-1,size(R1,2)));
    int_R3=(zeros(size(R3,1)*2-1,size(R3,2)));
    
    for i=1:size(int_R1,2)
        int_R1(:,i)=xcorr(R1(:,i),R1(:,ref));
        int_R3(:,i)=xcorr(R3(:,i),R3(:,ref));
    end
    
    figure(5);
    subplot(2,2,1)
    imagesc([1,Nr],[dt,dt*nt],R1);
    colorbar;
    title(num2str(ref));
    subplot(2,2,2)
    imagesc([1,Nr],[-dt*nt,dt*nt],int_R1);
    colorbar;
    subplot(2,2,3)
    imagesc([1,Nr],[dt,dt*nt],R3);
    colorbar;
    subplot(2,2,4)
    imagesc([1,Nr],[-dt*nt,dt*nt],int_R3);
    colorbar;
    shg;
end