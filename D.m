function A=D(nx,nz,direction)
%%
[i1,i2]=ind2sub([nx,nz],1:nx*nz);
A=sparse(nx*nz,nx*nz);
switch direction
    case 1
        A(sub2ind([nx*nz,nx*nz],1:nx*nz,sub2ind([nx,nz],i1,i2)))=-1;
        
        boundary=find(i1==nx);
        i1(boundary)=[];
        i2(boundary)=[];
        tt=1:nx*nz;
        tt(boundary)=[];
        A(sub2ind([nx*nz,nx*nz],tt,sub2ind([nx,nz],i1+1,i2)))=1;
    case -1
        A(sub2ind([nx*nz,nx*nz],1:nx*nz,sub2ind([nx,nz],i1,i2)))=1;
        
        boundary=find(i1==1);
        i1(boundary)=[];
        i2(boundary)=[];
        tt=1:nx*nz;
        tt(boundary)=[];
        A(sub2ind([nx*nz,nx*nz],tt,sub2ind([nx,nz],i1-1,i2)))=-1;
    case 3
        
        A(sub2ind([nx*nz,nx*nz],1:nx*nz,sub2ind([nx,nz],i1,i2)))=-1;
        
        boundary=find(i2==nz);
        i1(boundary)=[];
        i2(boundary)=[];
        tt=1:nx*nz;
        tt(boundary)=[];
        A(sub2ind([nx*nz,nx*nz],tt,sub2ind([nx,nz],i1,i2+1)))=1;
    case -3
        A(sub2ind([nx*nz,nx*nz],1:nx*nz,sub2ind([nx,nz],i1,i2)))=1;
        
        boundary=find(i2==1);
        i1(boundary)=[];
        i2(boundary)=[];
        tt=1:nx*nz;
        tt(boundary)=[];
        A(sub2ind([nx*nz,nx*nz],tt,sub2ind([nx,nz],i1,i2-1)))=-1;
end
end

