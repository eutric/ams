function [a,e,i,OM,om,th] = car2par(rr,vv,mu)

kk=[0 0 1]';
jj=[0 1 0]';
ii=[1 0 0]';
r=norm(rr);
v=norm(vv);
a=(2/r-v.^2/mu).^-1;
hh=cross(rr,vv);
h=norm(hh);
ee=cross(vv,hh)./mu-rr./r;
e=norm(ee);
i=acos(hh(3)./h);
N=cross(kk,hh)./norm(cross(kk,hh));
if N(2)>=0
    OM=acos(N(1));
else
    OM=2*pi-acos(N(1));
end
if ee(3)>=0
    om=acos(dot(N,ee)/e);
else
    om=2*pi-acos(dot(N,ee)/e);
end
vr=dot(vv,rr)/r;
if vr>=0
    th=acos(dot(rr,ee)/(r*e));
else
    th=2*pi-acos(dot(rr,ee)/(r*e));
end




end

