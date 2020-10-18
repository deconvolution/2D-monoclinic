function [C1,C2,C3,C7,C8,C12,C16,C19,C21,lambda2,mu2]=damp(C1,C2,C3,C7,C8,C12,C16,C19,C21,lambda2,mu2,l_damp,Q_interval,Q_cut)
% l_damp
Qs=[Q_cut,Q_cut+Q_interval:Q_interval:(Q_cut+Q_interval*(l_damp-1))];

if length(Qs)==1
    Qs=inf;
end
if length(Qs)==1
    Qs=inf;
end
%%
C1=C1;
C7=C7;
C8=C8;
C12=C12;
C16=C16;
C19=C19;
C21=C21;
C2=C2;
C3=C3;
lambda2=lambda2;
mu2=mu2;
[nx,ny,nz]=size(C1);
%% create cube for corners
cube=zeros(l_damp,l_damp,l_damp);
for i=1:l_damp
    cube(i,:,:)=Qs(l_damp-i+1);
    %cube(:,i,:)=Qs(l_damp-i+1);
    cube(:,:,i)=Qs(l_damp-i+1);
end
%% create plane
plane=zeros(l_damp,l_damp);
for i=1:l_damp
    plane(i,:)=Qs(l_damp-i+1);
    plane(:,i)=Qs(l_damp-i+1);
end
%% corner
cube2=ones(nx,ny,nz)*inf;
cube2(1:l_damp,1:l_damp,1:l_damp)=flip(flip(flip(cube,1),2),3);
cube2(end-l_damp+1:end,end-l_damp+1:end,end-l_damp+1:end)=cube;
cube2(end-l_damp+1:end,1:l_damp,1:l_damp)=flip(flip(cube,2),3);
cube2(1:l_damp,end-l_damp+1:end,1:l_damp)=flip(flip(cube,1),3);
cube2(1:l_damp,1:l_damp,end-l_damp+1:end)=flip(flip(cube,1),2);
cube2(1:l_damp,end-l_damp+1:end,end-l_damp+1:end)=flip(cube,1);
cube2(end-l_damp+1:end,1:l_damp,end-l_damp+1:end)=flip(cube,2);
cube2(end-l_damp+1:end,end-l_damp+1:end,1:l_damp)=flip(cube,3);
%%
cube2(1:l_damp,l_damp+1:end-l_damp,l_damp+1:end-l_damp)=repmat(reshape((Qs),[l_damp,1,1]),[1,ny-2*l_damp,nz-2*l_damp]);
cube2(end-l_damp+1:end,l_damp+1:end-l_damp,l_damp+1:end-l_damp)=repmat(reshape(flip(Qs),[l_damp,1,1]),[1,ny-2*l_damp,nz-2*l_damp]);
cube2(l_damp+1:end-l_damp,1:l_damp,l_damp+1:end-l_damp)=repmat(reshape((Qs),[1,l_damp,1]),[nx-2*l_damp,1,nz-2*l_damp]);
cube2(l_damp+1:end-l_damp,end-l_damp+1:end,l_damp+1:end-l_damp)=repmat(reshape(flip(Qs),[1,l_damp,1]),[nx-2*l_damp,1,nz-2*l_damp]);
cube2(l_damp+1:end-l_damp,l_damp+1:end-l_damp,1:l_damp)=repmat(reshape((Qs),[1,1,l_damp]),[nx-2*l_damp,ny-2*l_damp,1]);
cube2(l_damp+1:end-l_damp,l_damp+1:end-l_damp,end-l_damp+1:end)=repmat(reshape(flip(Qs),[1,1,l_damp]),[nx-2*l_damp,ny-2*l_damp,1]);
%
cube2(l_damp+1:end-l_damp,1:l_damp,1:l_damp)=repmat(reshape(flip(flip(plane,1),2),[1,l_damp,l_damp]),[nx-2*l_damp,1,1]);
cube2(l_damp+1:end-l_damp,end-l_damp+1:end,end-l_damp+1:end)=repmat(reshape(plane,[1,l_damp,l_damp]),[nx-2*l_damp,1,1]);
cube2(l_damp+1:end-l_damp,1:l_damp,end-l_damp+1:end)=repmat(reshape(flip(plane,1),[1,l_damp,l_damp]),[nx-2*l_damp,1,1]);
cube2(l_damp+1:end-l_damp,end-l_damp+1:end,1:l_damp)=repmat(reshape(flip(plane,2),[1,l_damp,l_damp]),[nx-2*l_damp,1,1]);

cube2(1:l_damp,l_damp+1:end-l_damp,1:l_damp)=repmat(reshape(flip(flip(plane,1),2),[l_damp,1,l_damp]),[1,ny-2*l_damp,1]);
cube2(end-l_damp+1:end,l_damp+1:end-l_damp,end-l_damp+1:end)=repmat(reshape(plane,[l_damp,1,l_damp]),[1,ny-2*l_damp,1]);
cube2(1:l_damp,l_damp+1:end-l_damp,end-l_damp+1:end)=repmat(reshape(flip(plane,1),[l_damp,1,l_damp]),[1,ny-2*l_damp,1]);
cube2(end-l_damp+1:end,l_damp+1:end-l_damp,1:l_damp)=repmat(reshape(flip(plane,2),[l_damp,1,l_damp]),[1,ny-2*l_damp,1]);

cube2(1:l_damp,1:l_damp,l_damp+1:end-l_damp)=repmat(reshape(flip(flip(plane,1),2),[l_damp,l_damp,1]),[1,1,nz-2*l_damp]);
cube2(end-l_damp+1:end,end-l_damp+1:end,l_damp+1:end-l_damp)=repmat(reshape(plane,[l_damp,l_damp,1]),[1,1,nz-2*l_damp]);
cube2(1:l_damp,end-l_damp+1:end,l_damp+1:end-l_damp)=repmat(reshape(flip(plane,1),[l_damp,l_damp,1]),[1,1,nz-2*l_damp]);
cube2(end-l_damp+1:end,1:l_damp,l_damp+1:end-l_damp)=repmat(reshape(flip(plane,2),[l_damp,l_damp,1]),[1,1,nz-2*l_damp]);
%%
cube3=zeros(nx,ny,nz);
cube3(cube2<inf)=1;
cube4=1-cube3;
C16=real(C16.*(1+1i/2./cube2).^2./(1+1/4./cube2.^2));
C19=real(C19.*(1+1i/2./cube2).^2./(1+1/4./cube2.^2));
C21=real(C21.*(1+1i/2./cube2).^2./(1+1/4./cube2.^2));
C1=real(C1.*(1+1i/2./cube2).^2./(1+1/4./cube2.^2));
C7=real(C7.*(1+1i/2./cube2).^2./(1+1/4./cube2.^2));
C12=real(C12.*(1+1i/2./cube2).^2./(1+1/4./cube2.^2));
C2=(C1-2*C16).*cube3+C2.*cube4;
C3=C2.*cube3+C3.*cube4;
C8=C2.*cube3+C8.*cube4;
mu2=imag(C16.*(1+1i./2./cube2).^2./(1+1/4./cube2.^2)).*cube3+mu2.*cube4;
lambda2=(imag(C12.*(1+1i/2./cube2).^2./(1+1/4./cube2.^2))-2*mu2).*cube3+lambda2.*cube4;
end