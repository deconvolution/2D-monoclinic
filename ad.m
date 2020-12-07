function u2=ad(u,n,direction)
switch direction
    case 1
        u2=(zeros(size(u,1)+n,size(u,2),'single'));
        u2(1:end-n,:)=u;
    case -1
        u2=(zeros(size(u,1)+n,size(u,2),'single'));
        u2(n+1:end,:)=u;
    case 3
        u2=(zeros(size(u,1),size(u,2)+n,'single'));
        u2(:,1:end-n)=u;
    case -3
        u2=(zeros(size(u,1),size(u,2)+n,'single'));
        u2(:,n+1:end)=u;
end