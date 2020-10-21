%% input
close all;
clear all;
dt=10^-3;
dx=8;
dz=8;
nt=1401;

nx=200;
nz=180;

r1=(2:10:(nx-1))*dx;
r3=(ones(size(r1))+20)*dz;

lp=20;
nPML=2;
R=.001;

rho=ones(nx,nz)*2000;

lambda=10^9;
mu=10^9;
C0=zeros(6,6);
C0(1,1)=lambda+2*mu;
C0(2,2)=C0(1,1);
C0(3,3)=C0(1,1);
C0(1,2)=lambda;
C0(1,3)=lambda;
C0(2,3)=lambda;
C0(4,4)=mu;
C0(5,5)=C0(4,4);
C0(6,6)=C0(4,4);

lambda2=zeros(nx,nz);
mu2=zeros(nx,nz);

C11=zeros(nx,nz);
C13=C11;
C15=C11;
C33=C11;
C35=C11;
C55=C11;

Eta11=C11;
Eta13=C11;
Eta15=C11;
Eta33=C11;
Eta35=C11;
Eta55=C11;

C11(:)=lambda+2*mu;
C13(:)=lambda;
%C15(:)=-.558*10^9;
C33(:)=lambda+2*mu;
%C35(:)=-.274*10^9;
C55(:)=mu;

scale=0;
Eta11=C11*scale;
Eta13=C13*scale;
Eta15=C15*scale;
Eta33=C33*scale;
Eta35=C35*scale;
Eta55=C55*scale;
% source
M=2.7;
s1=fix(nx/2)*dx;
s3=fix(nz/4)*dz;
freq=8;
singles=rickerWave(freq,dt,nt,M);
src1=zeros(nt,1);
src3=src1;
src1=1*[singles];
src3=1*[singles];
source_type='P';
plot_interval=200;
save_figure=1;
% figure path
p2=mfilename('fullpath');
path=[p2 '/'];
%%
[v1,v3,R1,R3]=monoclinic_2D_xz(dt,dx,dz,nt,nx,nz,...
    r1,r3,...
    s1,s3,src1,src3,source_type,...
    lp,nPML,R,...
    rho,C11,C13,C15,C33,C35,C55,...
    Eta11,Eta13,Eta15,Eta33,Eta35,Eta55,...
    plot_interval,...
    save_figure,path);