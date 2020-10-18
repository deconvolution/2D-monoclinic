%% input
close all;
clear all;

A=imread('NS.png');
A=im2double(A);
A=[repmat(A(1,:),[10,1]);
    A;
    repmat(A(end,:),[10,1])];
A=A(:,:);
solid=find(A==1);
fluid=find(A==0);
dt=10^-9;
dx=10^-5;
dz=10^-5;
nt=5000;

nx=size(A,1);
nz=size(A,2);

r1=(2:20:(nx-1))*dx;
r3=(ones(size(r1)))*(nz-11)*dz;

lp=10;
nPML=2;
R=.001;

rho=zeros(nx,nz);


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

mu=2.94*10^3*(1490)^2;
lambda=2.94*10^3*(2860)^2-2*mu;

C11(solid)=lambda;
C13(solid)=lambda;
C15(solid)=0;
C33(solid)=lambda+2*mu;
C35(solid)=0;
C55(solid)=mu;
rho(solid)=2.94*10^3;

C11(fluid)=1145*340^2;
C13(fluid)=1145*340^2;
C15(fluid)=0;
C33(fluid)=1145*340^2;
C35(fluid)=0;
C55(fluid)=0;
rho(fluid)=1145;

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
s3=11*dz;
sn=length(s1);
freq=5*10^6;
singles=rickerWave(freq,dt,nt,M);
src1=zeros(nt,1);
src3=src1;
src1=0*[singles];
src3=1*[singles];
source_type='D';
plot_interval=500;
save_figure=1;
path='./monoclinic_interface/';
%%
[v1,v3,R1,R3]=monoclinic_2D_xz(dt,dx,dz,nt,nx,nz,...
    r1,r3,...
    s1,s3,src1,src3,source_type,...
    lp,nPML,R,...
    rho,C11,C13,C15,C33,C35,C55,...
    Eta11,Eta13,Eta15,Eta33,Eta35,Eta55,...
    plot_interval,...
    save_figure,path);