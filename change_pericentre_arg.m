function [delta_v,th_i,th_f,th_best] = change_pericentre_arg(orbit,om_f,th_0)
a=orbit.a;
e=orbit.e;
om_i=orbit.om;
mu=orbit.mu;
D_om=om_f-om_i;
%teta di manovra in riferimento iniziale
th_1i=D_om/2;
th_2i=D_om/2+pi;
th_i=[th_1i;th_2i];
%teta di manovra in riferimento d arrivo
th_1f=2*pi-D_om/2;
th_2f=pi-D_om/2;
th_f=[th_1f;th_2f];
%calcolo p
p=a*(1-e^2);
%calcolo costo
delta_v=2*sqrt(mu/p)*e*sin(abs(D_om/2));
if D_om<=0
    disp('Delta omega nel trasferimento di pericentro è negativo')
end
%cerco teta migliore
if th_0>pi
    th_0 = pi-th_0;
end
if th_0>-D_om/2 && th_0<D_om/2
    th_best = D_om/2;
    disp("Il punto d'intersezione in cui effettuo la manovra è theta1")
else
    th_best = pi+D_om/2;
    disp("Il punto d'intersezione in cui effettuo la manovra è theta2")
end
end
