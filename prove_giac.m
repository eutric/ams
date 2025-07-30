%% Scenario 1
% Sequenze di manovre standard 
clear all
close all
clc

%parametri figura
dth=0.001;

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

[rr_start, vv_start] = par2car(O_start, th_start);

% Dati forniti dell'orbita d'arrivo
rr_end=[-12985.280000 3801.011400 8109.619300]';
vv_end=[-0.948600 -6.134000 1.356000]';

% Trovo dati mancanti dell'orbita d'arrivo 
[O_end,th_end] = car2par(rr_end,vv_end,mu);

%valutazione soluzione 1: standard
%1) piano orbita
[delta_v1, th_cp, orbit_cp] = changeOrbitalPlane(O_start, O_end); %la funzione calcola il costo solo in funzione del piano di arrivo, l'orbita che restituisce non è O_end
delta_t1 = TOF(O_start, th_start, th_cp);
%2)modifico anomalia pericentro
[delta_v2, th_i_cper, th_f_cper, th_best, orbit_chper] = change_pericentre_arg(orbit_cp, O_end.om-pi, th_cp);
delta_t2 = TOF(orbit_cp, th_cp, th_best);
%3)modifico orbita con ellittica bitangente 
[delta_v1_bt, delta_v2_bt, delta_t_bt, orbit_bt] = bitangentTransfer(orbit_chper, O_end, 'aa'); %testanto tutte 4 le possibilità conviene aa con la modifica del pericentro fatta prima 


DELTA_V_1 = delta_v1+delta_v2+delta_v1_bt+delta_v2_bt
DELTA_T_1 = delta_t1 + min(delta_t2) + delta_t_bt

sol_1=figure;
sol_1.Name="soluzione 1: cambio piano - cambio pericentro - trasferimento ellittico bitangente";
scatter3(0,0,0)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can
plotOrbit(O_start,0, 2*pi,dth,'b') % Orbita di partenza
plotOrbit(orbit_cp,0, 2*pi,dth,'r') % Orbita post cambio piano
plotOrbit(orbit_chper,th_best(2),pi,dth,'g')
plotOrbit(orbit_bt,0,pi,dth,'c')
plotOrbit(O_end,pi,0,dth,'m')
xlim([-1e5,1e5])
ylim([-1e5,1e5])
zlim([-1e5,1e5])
legend('attrattore','partenza','arrivo','sitema riferimento','orbita di inizio','orbita modificata di piano','orbita modificata anomalia pericentro','orbita biellittica ausiliaria','orbita finale')

