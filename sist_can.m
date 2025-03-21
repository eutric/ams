function [] = sist_can()
A=2*10^4*eye(3);
quiver3(zeros(1,3),zeros(1,3),zeros(1,3),A(1,:),A(2,:),A(3,:),'k');

end

