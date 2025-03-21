%% Scenario 1
clear 
close all
clc
% Parametri attrattore
mu=398600;

th0=0;
thf=2*pi;
dth=pi/100;

% Orbita 1
a_i=24400.0;
e_i=0.728300;
i_i=0.104700;
OM_i=1.514000;
om_i=3.107000;
th_i=1.665000;


% Numerose manovre


% Orbita d'Arrivo
rr=[-12985.280000 3801.011400 8109.619300]';
vv=[-0.948600 -6.134000 1.356000]';

[a_f,e_f,i_f,OM_f,om_f,th_f] = car2par(rr,vv,mu);

% Plot
plotOrbit(a_i,e_i,i_i,OM_i,om_i,th0,thf,dth,mu) % Plotto Orbita di Partenza
plotOrbit(a_f,e_f,i_f,OM_f,om_f,th0,thf,dth,mu) % Plotto Orbita d'Arrivo






