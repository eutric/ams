clear
close all
clc

O_start.a=24400.0;
O_start.e=0.728300;
O_start.i=0.104700;
O_start.OM=1.514000;
O_start.om=3.107000;
O_start.mu=398600;
th_start=1.665000;
[rr_i,vv_i] = par2car(O_start, th_start);
% Orbita d'Arrivo
rr=[-12985.280000 3801.011400 8109.619300]';
vv=[-0.948600 -6.134000 1.356000]';
% 
% [O_end,th_end] = car2par(rr,vv,mu);

e = 0.7;
E_1 = @(th)2*atan(sqrt((1-e)/(1+e)) .* tan(th/2));
ths = 0:0.01:2*pi;
figure
plot (ths, E_1(ths))
grid on

%% Ottimizzazione scenario 1
clear 
close all
clc
% Provo 1.Cambio Forma, 2. Cambio piano, 3. Cambio Pericentro, 4. attesa
mu = 398600;

O_start.a=24400.0;
O_start.e=0.728300;
O_start.i=0.104700;
O_start.OM=1.514000;
O_start.om=3.107000;
O_start.mu= mu;
th_start=1.665000;
[rr_i,vv_i] = par2car(O_start, th_start);
% Orbita d'Arrivo
rr=[-12985.280000 3801.011400 8109.619300]';
vv=[-0.948600 -6.134000 1.356000]';
% 
[O_end,th_end] = car2par(rr,vv,mu);

% Cambio Forma (provo tutti i casi, partire dal p dovrebbe essere più
% conveniente)
% type = 'pp';
% [delta_v1, delta_v2, delta_t, orbit_bt, th0, thf] = bitangentTransfer(O_start, O_end, type);
% costopp = delta_v1+delta_v2
% 
% type = 'pa';
% [delta_v1, delta_v2, delta_t, orbit_bt, th0, thf] = bitangentTransfer(O_start, O_end, type);
% costopa = delta_v1+delta_v2
% 
% type = 'ap';
% [delta_v1, delta_v2, delta_t, orbit_bt, th0, thf] = bitangentTransfer(O_start, O_end, type);
% costoap = delta_v1+delta_v2
% 
% type = 'aa';
% [delta_v1, delta_v2, delta_t, orbit_bt, th0, thf] = bitangentTransfer(O_start, O_end, type);
% costoaa = delta_v1+delta_v2
% Conviene pa

% 1.
type = 'pa';
[delta_v1, delta_v2, delta_t, orbit_bt, th0, thf, orbit_arr] = bitangentTransfer(O_start, O_end, type);
costopa = delta_v1+delta_v2

% 2.
[delta_v, om_f, theta, orbit_cp] = changeOrbitalPlane(orbit_arr, O_end);
costo_piano = delta_v

% 3.
[delta_v, th_i, th_f, th_best, orbit_chper] = change_pericentre_arg(orbit_cp, O_end.om, theta);
costo_pericentro = delta_v

costo_fin = costopa + costo_piano + costo_pericentro

