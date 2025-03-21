%%
clear all
close all
clc
mu=398600;
th0=0;
thf=2*pi;
dth=pi/100;

%orbita 1
a_i=24400.0;
e_i=0.728300;
i_i=0.104700;
OM_i=1.514000;
om_i=3.107000;
th_i=1.665000;

plotOrbit(a_i,e_i,i_i,OM_i,om_i,th0,thf,dth,mu)

%orbita finale

rr=[-12985.280000 3801.011400 8109.619300]';
vv=[-0.948600 -6.134000 1.356000]';

[a,e,i,OM,om,th] = car2par(rr,vv,mu);

plotOrbit(a,e,i,OM,om,th0,thf,dth,mu)






