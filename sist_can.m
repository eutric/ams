function [] = sist_can(sc)
if sc==1
A=2*10^4*eye(3);
quiver3(zeros(1,3),zeros(1,3),zeros(1,3),A(1,:),A(2,:),A(3,:),'k');
elseif sc==2
    A=2*10^8*eye(3);
quiver3(zeros(1,3),zeros(1,3),zeros(1,3),A(1,:),A(2,:),A(3,:),'k');
end
end

