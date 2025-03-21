function [orbit,th_] = car2par(rr,vv,mu)
% La funzione prende in ingresso i vettori posizione e velocità di un punto
% dell'orbita in coordinate cartesiane (sistema di riferimento Geocentrico
% inerziale) e restituisce i parametri caratterizzanti nel sistema di
% riferimento perifocale

% Input
% 
% Output

kk=[0 0 1]';
r=norm(rr);
v=norm(vv);
a_=(2/r-v.^2/mu).^-1;
hh=cross(rr,vv);
h=norm(hh);
ee=cross(vv,hh)./mu-rr./r;
e_=norm(ee);
i_=acos(hh(3)./h);
N=cross(kk,hh)./norm(cross(kk,hh));
if N(2)>=0
    OM_=acos(N(1));
else
    OM_=2*pi-acos(N(1));
end
if ee(3)>=0
    om_=acos(dot(N,ee)/e_);
else
    om_=2*pi-acos(dot(N,ee)/e_);
end
vr=dot(vv,rr)/r;
if vr>=0
    th_=acos(dot(rr,ee)/(r*e_));
else
    th_=2*pi-acos(dot(rr,ee)/(r*e_));
end
orbit.a=a_;
orbit.e=e_;
orbit.i=i_;
orbit.OM=OM_;
orbit.om=om_;
orbit.mu=mu;





end

