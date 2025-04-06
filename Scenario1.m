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

% Effettuo il trasferimento di orbita 
[delta_v1_bt, delta_v2_bt, delta_t_bt, orbit_bt] = bitangentTransfer(orbit_chper, O_end, 'pa');
plotOrbit(orbit_bt,th0,pi,dth,['-', 'c'])
% non penso abbiano molto senso i trasferimenti aa o pp per come sono messe
% le orbite 



% HO PROVATO A FAR GIRARE UN SATELLITE IN MANIERA APPROSSIMATIVA PERCHE' VOLEVO VEDERLO MUOVERSI 
% MA NON SI MUOVE COME DOVREBBE, VA AL CONTRARIO E NON SEGUE LEGGI FISICHE 

a = 24400.0;
e = 0.728300;
theta = linspace(0, 2*pi, 200); 
r = (a * (1 - e^2)) ./ (1 + e * cos(theta));
x = r .* cos(theta);
y = r .* sin(theta);
r2 = [x; y; zeros(length(x),1)'];


R_OM = [cos(O_start.OM) sin(O_start.OM) 0;
    -sin(O_start.OM) cos(O_start.OM) 0;
    0 0 1];
R_i = [1 0 0;
    0 cos(O_start.i) sin(O_start.i);
    0 -sin(O_start.i) cos(O_start.i)];
R_om = [cos(O_start.om) sin(O_start.om) 0;
    -sin(O_start.om) cos(O_start.om) 0;
    0 0 1];

T = R_om*R_i*R_OM;
rr2 = T' * r2; 


ss = plot3(rr2(1,:), rr2(2,:), rr2(3,:), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

for k = 1:length(theta)
    set(ss, 'XData', rr2(1,k), 'YData', rr2(2,k), 'ZData', rr2(3,k));

    pause(0.01); 
    drawnow; 
end

