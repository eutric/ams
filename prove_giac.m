%% Scenario 1
% Sequenze di manovre standard 
clear 
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


DELTA_V_1 = delta_v1 + delta_v2 + delta_v1_bt + delta_v2_bt
DELTA_T_1 = delta_t1 + min(delta_t2) + delta_t_bt

sol_1=figure;
sol_1.Name="STD: cambio piano - cambio pericentro - trasferimento ellittico bitangente";
scatter3(0,0,0)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can
plotOrbit(O_start,th_start,th_cp,dth,'b') % Orbita di partenza
plotOrbit(orbit_cp,th_cp,th_best(1),dth,'r') % Orbita post cambio piano
plotOrbit(orbit_chper,th_best(2),pi,dth,'g')
plotOrbit(orbit_bt,0,pi,dth,'c')
plotOrbit(O_end,pi,0,dth,'m')
xlim([-1e5,1e5])
ylim([-1e5,1e5])
zlim([-1e5,1e5])
legend('attrattore','partenza','arrivo','sitema riferimento','orbita di inizio','orbita modificata di piano','orbita modificata anomalia pericentro','orbita biellittica ausiliaria','orbita finale')


%% Scenario 1 - SOLUZIONE 02

clear 
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

%1)modifico orbita con bitangente 

%1) piano orbita
[delta_v1_bt, delta_v2_bt, delta_t, orbit_bt_temporanea, th0, thf, orbit_forma] = bitangentTransfer(O_start, O_end, 'aa');
%2) modifico piano
[delta_v2, th_cp, orbit_cp] = changeOrbitalPlane(orbit_forma, O_end);
delta_t2 = TOF(O_start, th_start, th_cp);
%3) cambio anomalia del pericentro
[delta_v3, th_i, th_f, th_best, orbit_chper] = change_pericentre_arg(orbit_cp, O_end.om, th_cp);
delta_t3 = TOF(orbit_cp, th_cp, th_best(1));
DELTA_V_2 = delta_v1_bt + delta_v2_bt + delta_v2 + delta_v3
DELTA_T_2 = delta_t + delta_t2 + delta_t3
%plot forma e poi attitudine
fig_2=figure;
fig_2.Name='bitangente - cambio piano - cambio pericentro';
scatter3(0,0,0)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can
plotOrbit(O_start, th_start, pi, dth, 'b') %           Start
plotOrbit(orbit_bt_temporanea, 0, pi, dth, 'c') %      Trasferimento
plotOrbit(orbit_forma, pi, th_cp, dth, 'g') %          O_end NON cambiata di piano
plotOrbit(orbit_cp, th_cp, th_best(1),dth,'k') %       Cambio Piano
plotOrbit(orbit_chper,th_best(2),th_end,dth,'m') % End
% plotOrbit(orbit_chper,th_best(2),pi,dth,'r')
legend('attrattore','partenza', 'arrivo', 'sistema ref','start','bitangente ausiliaria','orbita di forma','orbita di piano','orbita pericentro')

% Costo 5.0604

%% Scenario 1 - Soluzione 03
clear
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


% 1) Cambio forma
n = 1;
ra_t = n * O_start.a*(1+O_start.e); % n volte l'apocentro di O_start
[delta_v1_be, delta_v2_be, delta_v3_be, delta_t1, delta_t2, orbit_biel1, orbit_biel2] = biellipticTransfer(O_start,O_end, ra_t); % modifica a,e

% 2) Cambio piano
[delta_v2, th_cp, orbit_cp] = changeOrbitalPlane(orbit_biel2, O_end);

% 3) Cambio anomalia del pericentro
[delta_v3, th_i, th_f, th_best, orbit_chper] = change_pericentre_arg(orbit_cp, O_end.om, th_cp);

%3) finisco biellittico
orbit_biel2.i  = orbit_chper.i;
orbit_biel2.OM = orbit_chper.OM;
orbit_biel2.om = orbit_chper.om; % Riprendo la manovra, biel2 è cambiata

delta_t3 = TOF (orbit_biel1, pi, th_cp); % Attesa per cambiare piano
delta_t4 = TOF (orbit_cp, th_cp, th_best(1)); % Attesa per cambiare pericentro
delta_t5 = TOF (orbit_chper, th_best(2), pi); % Attesa per riprendere la manovra

% Poi si verifica delta_t2 lungo orbit_biel2

delta_t6 = TOF (O_end, 0, th_end);

DELTA_T3 = delta_t1 + delta_t3 + delta_t4 + delta_t5 + delta_t2 + delta_t6 % Penso manchi l'attesa iniziale

DELTA_V3 = delta_v1_be + delta_v2 + delta_v3 + delta_v2_be + delta_v3_be

fig_3=figure;
fig_3.Name='caso biellittico con modifica piano in apocentro';
scatter3(0,0,0)
hold on 
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can

plotOrbit(O_start,th_start, 0, dth, 'b') % Inizio
plotOrbit(orbit_biel1,0,th_cp,dth, 'm')  % Biel1 - fino al cambio piano
plotOrbit(orbit_cp, th_cp, th_best(1), dth, 'c')    % Biel1 - fino al cambio om
plotOrbit(orbit_chper, th_best(2), pi, dth, 'g')    % Biel1 torno all'apocentro
plotOrbit(orbit_biel2,pi,0,dth, 'k') % Biel2 fino a O_end
plotOrbit(O_end, 0, th_end, dth, 'r') % Fine

% plotOrbit(O_end,0,2*pi,dth, 'r')
legend ('terra','inizio', 'fine', 'SdR', 'start', 'prima biell','transfer al cambio piano','transfer al pericentro','seconda biell','end')

%COMMENTO: Manovre in Biel1 --> 8.3218
%          Manovre in Biel2 --> 4.9521

%% Scenario 1 - SOLUZIONE04

clear 
close all
% clc

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

%1)modifico orbita con bitangente 

%1) piano orbita

[dv_imp1, dv_imp2, delta_t, orbit_bt_temporanea, th0, thf, orbit_forma] = bitangentTransfer(O_start, O_end, 'aa');
%2) modifico piano
[dv_cp, th_cp, orbit_cp] = changeOrbitalPlane(orbit_bt_temporanea, O_end);
delta_t2 = TOF(O_start, th_start, th_cp);
%3) cambio anomalia del pericentro
[dv_cper, th_i, th_f, th_best, orbit_chper] = change_pericentre_arg(orbit_cp, O_end.om, th_cp);
delta_t3 = TOF(orbit_cp, th_cp, th_best(1));
DELTA_V_4 = abs(dv_imp1) + abs(dv_cp) + abs(dv_cper) + abs(dv_imp2)
dv_cp
dv_imp1
dv_imp2
dv_cper

DELTA_T_4 = delta_t + delta_t2 + delta_t3
%plot forma e poi attitudine
fig_2=figure;
fig_2.Name='bitangente - cambio piano - cambio pericentro';
scatter3(0,0,0)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can
plotOrbit(O_start, th_start, pi, dth, 'b') %           Start
plotOrbit(orbit_bt_temporanea, 0, th_cp, dth, 'c') %      Trasferimento
plotOrbit(orbit_cp, th_cp, th_best(1), dth, 'g') %          O_end NON cambiata di piano
plotOrbit(orbit_chper, th_best(2), pi,dth,'k') %       Cambio Piano
plotOrbit(O_end, pi, th_end,dth,'m') % End
% plotOrbit(orbit_chper,th_best(2),pi,dth,'r')
legend('attrattore','partenza', 'arrivo', 'sistema ref','start','bitangente ausiliaria','orbita di forma','orbita di piano','orbita pericentro')

% Costo 4.2537
figure
plotOrbit_plane(orbit,th0,thf,dth,linestyle)

%% Scenario 1 - SOLUZIONE 05
clear
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

%1)modifico orbita con bitangente 

%1) piano orbita
[dv_cp, th_cp, orbit_cp] = changeOrbitalPlane(O_start, O_end);
[dv_imp1, dv_imp2, delta_t, orbit_bt_temporanea, th0, thf, orbit_forma] = bitangentTransfer(orbit_cp, O_end, 'aa');
%2) modifico piano

delta_t2 = TOF(O_start, th_start, th_cp);
%3) cambio anomalia del pericentro
[dv_cper, th_i, th_f, th_best, orbit_chper] = change_pericentre_arg(orbit_bt_temporanea, O_end.om, th_cp);
delta_t3 = TOF(orbit_cp, th_cp, th_best(1));
DELTA_V_4 = abs(dv_imp1) + abs(dv_cp) + abs(dv_cper) + abs(dv_imp2)
dv_cp
dv_imp1
dv_imp2
dv_cper

DELTA_T_4 = delta_t + delta_t2 + delta_t3
%plot forma e poi attitudine
fig_2=figure;
fig_2.Name='bitangente - cambio piano - cambio pericentro';
scatter3(0,0,0)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can
plotOrbit(O_start, th_start, th_cp, dth, 'b') %           Start
plotOrbit(orbit_cp, th_cp, pi, dth, 'c') %      Trasferimento
plotOrbit(orbit_bt_temporanea, 0, th_best(1), dth, 'g') %          O_end NON cambiata di piano
plotOrbit(orbit_chper, th_best(2), pi,dth,'k') %       Cambio Piano
plotOrbit(O_end, pi, th_end,dth,'m') % End
% plotOrbit(orbit_chper,th_best(2),pi,dth,'r')
legend('attrattore','partenza', 'arrivo', 'sistema ref','start','bitangente ausiliaria','orbita di forma','orbita di piano','orbita pericentro')

% Costo 4.2537
figure
subplot(1,2,1)

scatter(0,0)
hold on
xline(0)
yline(0)
plotOrbit_plane(O_start, th_start, th_cp, dth, 'b')
legend('attrattore','x', 'y','vado al cambio piano')
grid on
subplot(1,2,2)
scatter(0,0)
hold on
xline(0)
yline(0)
% plotOrbit_plane(orbit_cp, th_start, th_cp, dth, 'b')
plotOrbit_plane(orbit_cp, th_cp, pi, dth, 'c')
plotOrbit_plane(orbit_bt_temporanea, 0, th_best(1), dth, 'g')
plotOrbit_plane(orbit_chper, th_best(2), pi, dth, 'k')
plotOrbit_plane(O_end, pi, th_end, dth, 'm')
plotOrbit_plane(O_end, 0, 2*pi, dth, '--')
grid on
legend('attrattore','x', 'y','vado al primo impulso','trasferimento','vado al secondo impulso','fine')

%% ricerca del miglior cambio piano
clear
close all
clc

% hx*x+hy*y+hz*z = 0, z = (-hx*x-hy*y)/hz
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

z = @(h,x,y) (-h(1)*x -h(2)*y)/h(3);
x = linespace (-10^4,10^4, 100);
y = linespace (-10^4,10^4, 100);
[X, Y] = meshgrid (x,y);
Z = z(h, X, Y);

surf (X, Y, Z)

%% Sistemo change_pericenter
clear
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

%1)modifico orbita con bitangente 

%1) piano orbita
[dv_cp, th_cp, orbit_cp] = changeOrbitalPlane(O_start, O_end);
[dv_imp1, dv_imp2, delta_t, orbit_bt_temporanea, th0, thf, orbit_forma] = bitangentTransfer(orbit_cp, O_end, 'aa');
%2) modifico piano

delta_t2 = TOF(O_start, th_start, th_cp);
%3) cambio anomalia del pericentro
[dv_cper, th_i, th_f, th_best, orbit_chper] = change_pericentre_arg(orbit_bt_temporanea, O_end.om, th_cp);

delta_t3 = TOF(orbit_cp, th_cp, th_best(1));
DELTA_V_4 = abs(dv_imp1) + abs(dv_cp) + abs(dv_cper) + abs(dv_imp2)
% Costo 4.2537
figure
subplot(1,2,1)

scatter(0,0)
hold on
xline(0)
yline(0)
plotOrbit_plane(O_start, th_start, th_cp, dth, 'b')
legend('attrattore','x', 'y','vado al cambio piano')
grid on
subplot(1,2,2)
scatter(0,0)
hold on
xline(0)
yline(0)
% plotOrbit_plane(orbit_cp, th_start, th_cp, dth, 'b')
plotOrbit_plane(orbit_cp, th_cp, pi, dth, 'c')
plotOrbit_plane(orbit_bt_temporanea, 0, th_best(1), dth, 'g')
plotOrbit_plane(orbit_chper, th_best(2), pi, dth, 'k')
plotOrbit_plane(O_end, pi, th_end, dth, 'm')
plotOrbit_plane(O_end, 0, 2*pi, dth, '--')
grid on
legend('attrattore','x', 'y','vado al primo impulso','trasferimento','vado al secondo impulso','fine')

%% Scenario 2
clear
close all
clc

mu = 1.32712440042e20 * 0.001^3; % km^3/s^2
dth = 0.001;

% Terra
O_start.a = 1.4946e8;  % [km]
O_start.e = 0.016;     % [ ]
O_start.i = 9.1920e-5; % [rad]
O_start.OM = 2.7847;   % [rad]
O_start.om = 5.2643;   % [rad]
O_start.mu = mu;       % [km^3/s^2]

% Asteroide 163899 (2003 SD220)
O_end.a = 0.827903 *1.496e8; % [km]
O_end.e = 0.209487;          % [ ]
O_end.i = 8.55 *pi/180;      % [rad]
O_end.OM = 273.63 *pi/180;   % [rad]
O_end.om = 327.03 *pi/180;   % [rad]
O_end.mu = mu;       
