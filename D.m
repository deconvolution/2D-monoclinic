function u2=D(u,direction)
%%
switch direction
    case 1
        u2=-1*ad(u(1:end-1,:),1,-1) ...
            +1*u;
    case -1
        u2=-1*u ...
            +1*ad(u(2:end,:),1,1);
    case 3
        u2=-1*ad(u(:,1:end-1),1,-3) ...
            +1*u;
    case -3
        u2=-1*u ...
            +1*ad(u(:,2:end),1,3);
end
end

