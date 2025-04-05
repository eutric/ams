%% Scenario 1
% Sequenze di manovre standard 
clear 
close all
clc


% Parametro attrattore
mu=398600;
% Dati forniti dell'orbita di partenza
O_start.a = 24400.0;
O_start.e = 0.728300;
O_start.i = 0.104700;
O_start.OM = 1.514000;
O_start.om = 3.107000;
O_start.mu = mu;
th_start = 1.665000;


% Dati forniti dell'orbita d'arrivo
rr=[-12985.280000 3801.011400 8109.619300]';
vv=[-0.948600 -6.134000 1.356000]';

% Trovo dati mancanti dell'orbita d'arrivo 
[O_end,th_end] = car2par(rr,vv,mu);

% Plotto le due orbite
th0=0;
thf=2*pi;
dth=pi/100;

plotOrbit(O_start,th0,thf,dth, 'b') % Partenza
hold on
plotOrbit(O_end,th0,thf,dth, 'k') % Arrivo

% Plotto il pianeta Terra
[X,Y,Z]=sphere(1000);
k=0.5*10^4;
surf(k*X,k*Y,k*Z,"LineStyle","none","FaceColor","texturemap",CData=flip((imread('earth.jpg'))));
xlim([-0.6*10^5;0.6*10^5])
ylim([-0.6*10^5;0.6*10^5])
zlim([-0.6*10^5;0.6*10^5])

% Modifico il piano dell'orbita iniziale 
[delta_v1, om_cp, th_cp, orbit_cp] = changeOrbitalPlane(O_start, O_end);
plotOrbit(orbit_cp,th0,thf,dth,['--', 'r'])

% Modifico anomalia del pericentro
[delta_v2, th_i, th_f, th_best, orbit_chper] = change_pericentre_arg(orbit_cp, O_end.om, th_cp);
plotOrbit(orbit_chper,th0,thf,dth,['--', 'g'])


[delta_v1_bt, delta_v2_bt, delta_t_bt] = bitangentTransfer(orbit_chper, O_end, 'pa');


