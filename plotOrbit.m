function [] = plotOrbit(orbit,th0,thf,dth)
a=orbit.a;
e=orbit.e;
i=orbit.i;
om=orbit.om;
OM=orbit.OM;
mu=orbit.mu;

teta=th0:dth:thf;
p=a*(1-e^2);
r=@(t)p./(1+e.*cos(t));
rv=r(teta);
rx=rv.*cos(teta);
ry=rv.*sin(teta);
r2=[rx;ry;zeros(length(rx),1)'];
R_OM=[cos(OM) sin(OM) 0;
    -sin(OM) cos(OM) 0;
    0 0 1];
R_i=[1 0 0;
    0 cos(i) sin(i);
    0 -sin(i) cos(i)];
R_om=[cos(om) sin(om) 0;
    -sin(om) cos(om) 0;
    0 0 1];
T=R_om*R_i*R_OM;
rr=T'*r2;
plot3(rr(1,:),rr(2,:),rr(3,:))


grid on
end

