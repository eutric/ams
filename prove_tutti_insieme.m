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
O_start.a=24400.0;
O_start.e=0.728300;
O_start.i=0.104700;
O_start.OM=1.514000;
O_start.om=3.107000;
O_start.mu=mu;
th_start=1.665000;


% Numerose manovre


% Orbita d'Arrivo
rr=[-12985.280000 3801.011400 8109.619300]';
vv=[-0.948600 -6.134000 1.356000]';

[O_end,th] = car2par(rr,vv,mu);

% Plot
plotOrbit(O_start,th0,thf,dth)
hold on% Plotto Orbita di Partenza
plotOrbit(O_end,th0,thf,dth) % Plotto Orbita d'Arrivo

%% modifica orbita anomalia del pericentro
clear all
close all
clc
%dati attrattore
mu=398600;
%dati orbita inizio in forma parametrica
O_start.a=24400.0;
O_start.e=0.728300;
O_start.i=0.104700;
O_start.OM=1.514000;
O_start.om=3.107000;
O_start.mu=mu;
th_start=1.665000;


%dati orbita arrivo in forma parametrica
rr=[-12985.280000 3801.011400 8109.619300]';
vv=[-0.948600 -6.134000 1.356000]';
[O_end,th] = car2par(rr,vv,mu);

%primo orbita di modifica w

[delta_v,th_i,th_f,th_best] = change_pericentre_arg(O_start,O_end.om,th_start);

%generazione orbita con w uguale
O_w=O_start;
O_w.om=O_end.om;
%plot orbite
th0=0;
thf=2*pi;
dth=pi/100;
plotOrbit(O_start,th0,thf,dth) % Plotto Orbita di Partenza
hold on
plotOrbit(O_w,th0,thf,dth)
legend('orbita1','orbita2')
sist_can
[X,Y,Z]=sphere(1000);
k=0.5*10^4;
surf(k*X,k*Y,k*Z,"LineStyle","none","FaceColor","texturemap",CData=flip((imread('pepera.jpg'))));
xlim([-10^5;10^5])
ylim([-10^5;10^5])
zlim([-10^5;10^5])






