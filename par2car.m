function [rr,vv] = par2car(a,e,i,OM,om,th,mu)

R_OM=[cos(OM) sin(OM) 0;
    -sin(OM) cos(OM) 0;
    0 0 1];
R_i=[1 0 0;
    0 cos(i) sin(i);
    0 -sin(i) cos(i)];
R_om=[cos(om) sin(om) 0;
    -sin(om) cos(om) 0;
    0 0 1];
p=a*(1-e^2);
r=p/(1+e*cos(th));
rr2=r*[cos(th); sin(th);0];%2 xke in 2D
vv2=sqrt(mu/p)*[-sin(th);e+cos(th);0];
T=R_om*R_i*R_OM;
rr=T'*rr2;
vv=T'*vv2;
end

