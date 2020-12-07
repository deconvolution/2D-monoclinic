function [v1,v3,R1,R3]=monoclinic_2D_xz(dt,dx,dz,nt,nx,nz,...
    r1,r3,...
    s1,s3,src1,src3,source_type,...
    lp,nPML,R,...
    rho,C11,C13,C15,C33,C35,C55,...
    Eta11,Eta13,Eta15,Eta33,Eta35,Eta55,...
    plot_interval,plot_source,...
    save_figure,path)
%% create folder for figures
if ~exist(path,'dir')
    mkdir(path)
end
n_picture=1;
%% PML
beta0=(ones(nx,nz,'single').*1000*(nPML+1)*log(1/R)/2/lp/dx);
beta1=(zeros(nx,nz,'single'));
beta3=beta1;
tt=(1:lp)/lp;
tt2=repmat(reshape(tt,[lp,1]),[1,nz]);
plane_grad1=zeros(nx,nz);
plane_grad3=plane_grad1;

plane_grad1(2:lp+1,:)=flip(tt2,1);
plane_grad1(nx-lp:end-1,:)=tt2;
plane_grad1(1,:)=plane_grad1(2,:);
plane_grad1(end,:)=plane_grad1(end-1,:);

tt2=repmat(reshape(tt,[1,lp]),[nx,1]);
plane_grad3(:,2:lp+1)=flip(tt2,2);
plane_grad3(:,nz-lp:end-1)=tt2;
plane_grad3(:,1)=plane_grad3(:,2);
plane_grad3(:,end)=plane_grad3(:,end-1);

beta1=beta0.*plane_grad1.^nPML;
beta3=beta0.*plane_grad3.^nPML;
%% source
ind_sor=sub2ind([nx,nz],s1,s3);
%% model configuration
R1=(zeros(nt,length(r1),'single'));
R3=(zeros(nt,length(r3),'single'));
ind_rec=sub2ind([nx,nz,3],r1,r3,ones(1,length(r1))*3);
%% monoclinic 2D solver xz plane (symmetric plane)
v1=(zeros(nx,nz,3,'single'));
v3=v1;
lim2=(zeros(4,2,'single'));

sigmas11=(zeros(nx,nz,'single'));
sigmas13=sigmas11;
sigmas33=sigmas11;
p=sigmas11;
ts=sigmas11;
ts2=sigmas11;

E=sigmas11;
F=sigmas11;
G=sigmas11;
H=sigmas11;
I=sigmas11;
J=sigmas11;
K=sigmas11;
L=sigmas11;
M=sigmas11;
N=sigmas11;
O=sigmas11;
P=sigmas11;
Q=sigmas11;
R=sigmas11;
S=sigmas11;

l=2;
ts(ind_sor)=src1(1,:);
ts2(ind_sor)=src3(1,:);

switch source_type
    case 'D'
        v1(:,:,3)=v1(:,:,3)+.5./rho.*ts;
        v3(:,:,3)=v3(:,:,3)+.5./rho.*ts2;
    case 'P'
        v1(:,:,3)=v1(:,:,3)+.5./rho.*D(ts,-1)/dx;
        v3(:,:,3)=v3(:,:,3)+.5./rho.*D(ts2,-3)/dz;
end
%%
tic;
for l=2:nt-1
    for l2=1:2
        v1(:,:,l2)=v1(:,:,l2+1);
        v3(:,:,l2)=v3(:,:,l2+1);
    end
    %% compute sigma
    v1_x=D(v1(:,:,2),1)/dx;
    v3_x=D(v3(:,:,2),-1)/dx;
    v1_z=D(v1(:,:,2),-3)/dz;
    v3_z=D(v3(:,:,2),3)/dz;
    
    v1_x2=D(v1(:,:,1),1)/dx;
    v3_x2=D(v3(:,:,1),-1)/dx;
    v1_z2=D(v1(:,:,1),-3)/dz;
    v3_z2=D(v3(:,:,1),3)/dz;
    
    E=dt*v1_x+E;
    F=dt*v1_z+F;
    G=dt*v3_x+G;
    H=dt*v3_z+H;
    
    I=dt*sigmas11+I;
    J=dt*sigmas33+J;
    K=dt*sigmas13+K;
    L=dt*p+L;
    
    sigmas11=dt*.5*((C11-C13).*v1_x+beta3.*(C11-C13).*E...
        +(C15-C35).*v3_x+beta3.*(C15-C35).*G...
        +(C13-C33).*v3_z+beta1.*(C13-C33).*H...
        +(C15-C35).*v1_z+beta1.*(C15-C35).*F...
        +(Eta11-Eta13)/dt.*(v1_x-v1_x2)+(Eta11-Eta13).*beta3.*v1_x...
        +(Eta15-Eta35)/dt.*(v3_x-v3_x2)+(Eta15-Eta35).*beta3.*v3_x...
        +(Eta13-Eta33)/dt.*(v3_z-v3_z2)+(Eta13-Eta33).*beta1.*v3_z...
        +(Eta15-Eta35)/dt.*(v1_z-v1_z2)+(Eta15-Eta35).*beta1.*v1_z)...
        +sigmas11...
        -dt*(beta1+beta3).*sigmas11...
        -dt*beta1.*beta3.*I;
    
    sigmas33=dt*.5*((C13-C11).*v1_x+beta3.*(C13-C11).*E...
        +(C35-C15).*v3_x+beta3.*(C35-C15).*G...
        +(C33-C13).*v3_z+beta1.*(C33-C13).*H...
        +(C35-C15).*v1_z+beta1.*(C35-C15).*F...
        +(Eta13-Eta11)/dt.*(v1_x-v1_x2)+(Eta13-Eta11).*beta3.*v1_x...
        +(Eta35-Eta15)/dt.*(v3_x-v3_x2)+(Eta35-Eta15).*beta3.*v3_x...
        +(Eta33-Eta13)/dt.*(v3_z-v3_z2)+(Eta33-Eta13).*beta1.*v3_z...
        +(Eta35-Eta15)/dt.*(v1_z-v1_z2)+(Eta35-Eta15).*beta1.*v1_z)...
        +sigmas33...
        -dt*(beta1+beta3).*sigmas33...
        -dt*beta1.*beta3.*J;
    
    sigmas13=dt*(C15.*v1_x+beta3.*C15.*E...
        +C35.*v3_z+beta1.*C35.*H...
        +C55.*(v1_z+v3_x)+C55.*(beta1.*F+beta3.*G)...
        +Eta15/dt.*(v1_x-v1_x2)+Eta15.*beta3.*v1_x...
        +Eta35/dt.*(v3_z-v3_z2)+Eta35.*beta1.*v3_z...
        +Eta55/dt.*(v1_z-v1_z2+v3_x-v3_x2)+Eta55.*(beta1.*v1_z+beta3.*v3_x))...
        +sigmas13...
        -dt*(beta1+beta3).*sigmas13...
        -dt*beta1.*beta3.*K;
    %% compute p
    p=dt*(-.5)*((C11+C13).*v1_x+beta3.*(C11+C13).*E...
        +(C13+C33).*v3_z+beta1.*(C13+C33).*H...
        +(C15+C35).*(v1_z+v3_x)+(C15+C35).*(beta1.*F+beta3.*G)...
        +(Eta11+Eta13)/dt.*(v1_x-v1_x2)+(Eta11+Eta13).*beta3.*v1_x...
        +(Eta13+Eta33)/dt.*(v3_z-v3_z2)+(Eta13+Eta33).*beta1.*v3_z...
        +(Eta15+Eta35)/dt.*(v1_z-v1_z2+v3_x-v3_x2)+(Eta15+Eta35).*(beta1.*v1_z+beta3.*v3_x))...
        +p...
        -dt*(beta1+beta3).*p...
        -dt*beta1.*beta3.*L;
    %% compute v
    M=dt*D(sigmas11-p,-1)/dx+M;
    N=dt*D(sigmas13,3)/dz+N;
    P=dt*D(sigmas13,1)/dx+P;
    Q=dt*D(sigmas33-p,-3)/dz+Q;
    R=dt*v1(:,:,2)+R;
    S=dt*v3(:,:,2)+S;
    v1(:,:,3)=dt./rho.*(D(sigmas11-p,-1)/dx+beta3.*M...
        +D(sigmas13,3)/dz+beta1.*N)...
        +v1(:,:,3)...
        -dt*(beta1+beta3).*v1(:,:,2)...
        -dt*beta1.*beta3.*R;
    v3(:,:,3)=dt./rho.*(D(sigmas13,1)/dx+beta3.*P...
        +D(sigmas33-p,-3)/dz+beta1.*Q)...
        +v3(:,:,3)...
        -dt*(beta1+beta3).*v3(:,:,2)...
        -dt*beta1.*beta3.*S;
    
    ts(ind_sor)=src1(l,:);
    ts2(ind_sor)=src3(l,:);
    
    switch source_type
        case 'D'
            v1(:,:,3)=v1(:,:,3)+1./rho.*ts;
            v3(:,:,3)=v3(:,:,3)+1./rho.*ts2;
        case 'P'
            v1(:,:,3)=v1(:,:,3)+1./rho.*D(ts,-1)/dx;
            v3(:,:,3)=v3(:,:,3)+1./rho.*D(ts2,-3)/dz;
    end
    %% fixed boundary condition
    v1(1,:,3)=0;
    v1(end,:,3)=0;
    v1(:,1,3)=0;
    v1(:,end,3)=0;
    
    v3(1,:,3)=0;
    v3(end,:,3)=0;
    v3(:,1,3)=0;
    v3(:,end,3)=0;
    %% assign recordings
    R1(l+1,:)=v1(ind_rec);
    R3(l+1,:)=v3(ind_rec);
    %% plot
    if mod(l,plot_interval)==0 || l==nt-1
        lim2(1,1)=min(min(v3,[],'all'),lim2(1,1));
        lim2(1,2)=max(max(v3,[],'all'),lim2(1,2));
        lim2(2,1)=min(min(sigmas33,[],'all'),lim2(2,1));
        lim2(2,2)=max(max(sigmas33,[],'all'),lim2(2,2));
        lim2(3,1)=min(min(sigmas13,[],'all'),lim2(3,1));
        lim2(3,2)=max(max(sigmas13,[],'all'),lim2(3,2));
        lim2(4,1)=min(min(p,[],'all'),lim2(4,1));
        lim2(4,2)=max(max(p,[],'all'),lim2(4,2));
        
        
        hfig=figure('Visible','off');
        set(gcf,'position',[80,80,1300,600]);
        subplot(2,3,1)
        imagesc([1,nx]*dx,[1,nz]*dz,v3(:,:,3)',.02*[lim2(1,1),lim2(1,2)+10^-12]);
        colorbar;
        xlabel({['x [m]']});
        ylabel({['z [m]']});
        title({['t=' num2str(l*dt) 's'],['v_3 [m/s]']});
        xlabel('x [m]');
        ylabel('z [m]');
        colorbar;
        hold on;
        if plot_source==1
            for i=1:length(s1)
                ax2=plot(s1(i)*dx,s3(i)*dz,'v','color',[1,0,0]);
                hold on;
            end
        end
        for i=1:length(r1)
            ax4=plot(r1(i)*dx,r3(i)*dz,'^','color',[0,1,1]);
            hold on;
        end
        ax3=plot([lp+1,lp+1]*dx,[lp+1,nz-lp-1]*dz,'color','blue');
        hold on;
        ax3=plot([nx-lp-1,nx-lp-1]*dx,[lp+1,nz-lp-1]*dz,'color','blue');
        hold on;
        hold on;
        ax3=plot([lp+1,nx-lp-1]*dx,[nz-lp-1,nz-lp-1]*dz,'color','blue');
        axis on;
        ax3=plot([lp+1,nx-lp-1]*dx,[lp+1,lp+1]*dz,'color','blue');
        axis on;
        
        subplot(2,3,2)
        imagesc([1,nx]*dx,[1,nz]*dz,sigmas33',.02*[lim2(2,1),lim2(2,2)+10^-12]);
        xlabel('x [m]');
        ylabel('z [m]');
        title('\sigma_{s33} [Pa]');
        colorbar;
        
        subplot(2,3,4)
        imagesc([1,nx]*dx,[1,nz]*dz,p',.02*[lim2(4,1),lim2(4,2)+10^-12]);
        xlabel('x [m]');
        ylabel('z [m]');
        title('p [Pa]');
        colorbar;
        
        subplot(2,3,3)
        imagesc([1,nx]*dx,[1,nz]*dz,sigmas13',.02*[lim2(3,1),lim2(3,2)+10^-12]);
        xlabel('x [m]');
        ylabel('z [m]');
        title('\sigma_{s13} [Pa]');
        colorbar;
        
        subplot(2,3,5)
        imagesc([1,length(r1)],[1,(l+1)]*dt,R3(1:l+1,:));
        colorbar;
        xlabel('Nr');
        ylabel('t [s]');
        title('R3 [m/s]');
        ylim([1,nt]*dt);
        
        subplot(2,3,6)
        imagesc([1,nx]*dx,[1,nz]*dz,C33');
        xlabel('x [m]');
        ylabel('z [m]');
        title('C33 [Pa]');
        colorbar;
        hold on;
        if plot_source==1
            for i=1:length(s1)
                ax2=plot(s1(i)*dx,s3(i)*dz,'v','color',[1,0,0]);
                hold on;
            end
        end
        for i=1:length(r1)
            ax4=plot(r1(i)*dx,r3(i)*dz,'^','color',[0,1,1]);
            hold on;
        end
        ax3=plot([lp+1,lp+1]*dx,[lp+1,nz-lp-1]*dz,'color','blue');
        hold on;
        ax3=plot([nx-lp-1,nx-lp-1]*dx,[lp+1,nz-lp-1]*dz,'color','blue');
        hold on;
        hold on;
        ax3=plot([lp+1,nx-lp-1]*dx,[nz-lp-1,nz-lp-1]*dz,'color','blue');
        axis on;
        ax3=plot([lp+1,nx-lp-1]*dx,[lp+1,lp+1]*dz,'color','blue');
        axis on;
        if plot_source==1
            legend([ax2,ax3,ax4],...
                'source','PML boundary','receiver',...
                'Location',[0.5,0.02,0.005,0.002],'orientation','horizontal');
        else
            legend([ax3,ax4],...
                'PML boundary','receiver',...
                'Location',[0.5,0.02,0.005,0.002],'orientation','horizontal');
        end
        if save_figure==1
            print(gcf,[path num2str(n_picture) '.png'],'-dpng','-r200');
            n_picture=n_picture+1;
        end
    end
    fprintf('\n time step=%d/%d',l+1,nt);
    fprintf('\n    epalsed time=%.2fs',toc);
    fprintf('\n    n_picture=%d',n_picture);
    d=clock;
    fprintf('\n    current time=%d %d %d %d %d %.0d',d(1),d(2),d(3),d(4),d(5),d(6));
end
%%
