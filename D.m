function u2=D(u,direction)
%%
switch direction
    case 1
        u2=[u;
            zeros(1,size(u,2))];
        u2=-u2(1:end-1,:)+u2(2:end,:);
    case -1
        u2=[zeros(1,size(u,2));
            u];
        u2=-1*u2(1:end-1,:)+u2(2:end,:);
    case 3
        u2=[u,zeros(size(u,1),1)];
        u2=-1*u2(:,1:end-1)+u2(:,2:end);
    case -3
        u2=[zeros(size(u,1),1),u];
        u2=-1*u2(:,1:end-1)+u2(:,2:end);
end
end

