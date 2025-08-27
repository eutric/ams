%% ricerca del min dv
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
O_end.mu = mu;               % [km^3/s^2]

fun_ob_dv = @(x) O_tfun(O_start, O_end, x(1), x(2), x(3), [0,0]).cost;

% Approccio a tappeto
n = 10; % Dimensioni griglia
a = [3.4, 3.4, 3.4]; % Si possono cambiare ad ogni lancio per infittire dove serve
b = [3.8, 3.8, 3.8]; 
[O_best_tappeto_dv] = tappeto(O_start, O_end, fun_ob_dv, n, a, b);

% ga
x_ga_dv = ga(fun_ob_dv,3,[],[],[],[],a,b);
O_best_ga_dv = O_tfun(O_start, O_end,x_ga_dv(1),x_ga_dv(2),x_ga_dv(3),[0,0]);

% fmincon
x0 = (a-b)./2;
x_fmc_dv = fmincon(fun_ob_dv,x0,[],[],[],[],a,b);
O_best_fmc_dv = O_tfun(O_start, O_end,x_fmc_dv(1),x_fmc_dv(2),x_fmc_dv(3),[0,0]);

fun_ob_dt = @(x) O_tfun(O_start, O_end, x(1), x(2), x(3), [0,0]).tempo;

% Tappeto
n = 20; % Dimensioni griglia
a = [5.9, 3.5, 1.7]; % Si possono cambiare ad ogni lancio per infittire dove serve
b = [6.1, 3.7, 1.8];
[O_best_tappeto_dt] = tappeto(O_start, O_end, fun_ob_dt, n, a, b);

% ga
x_ga_dt = ga(fun_ob_dt, 3, [], [], [], [], a, b);
O_best_ga_dt = O_tfun(O_start, O_end,x_ga_dt(1),x_ga_dt(2),x_ga_dt(3), [0,0]);

% fmincon
x0 = (a-b)./2;
x_fmc_dt = fmincon(fun_ob_dt, x0, [], [], [], [], a, b);
O_best_fmc_dt = O_tfun(O_start, O_end, x_fmc_dt(1), x_fmc_dt(2), x_fmc_dt(3), [0,0]);

O_best_dv = O_best_fmc_dv
O_best_dt = O_best_fmc_dt

% Punti di trasferimento
[rr1_mindeltav, vv1_mindeltav] = par2car (O_best_dv, O_best_dv.th_t(1));
[rr2_mindeltav, vv2_mindeltav] = par2car (O_best_dv, O_best_dv.th_t(2));

[rr1_mindeltat, vv1_mindeltat] = par2car (O_best_dt, O_best_dt.th_t(1));
[rr2_mindeltat, vv2_mindeltat] = par2car (O_best_dt, O_best_dt.th_t(2));

% plots

% orbite a confronto
figure
scatter3 (0,0,0, 100, [1, 0.5, 0], 'filled')
hold on
plotOrbit (O_start, 0, 2*pi, dth, 'b');
plotOrbit (O_end, 0, 2*pi, dth, 'r');
legend ('SOLE', 'Orbita Terrestre', 'Orbita Asteroide')

% min dv
figure_dv = figure;
scatter3 (0,0,0,100,"yellow", 'filled')
hold on
plotOrbit(O_start, 0, 2*pi, dth, 'b--');
plotOrbit(O_end, 0, 2*pi, dth, 'r--');
plotOrbit(O_best_dv, O_best_dv.th_t(1), O_best_dv.th_t(2), dth, 'g');
scatter3 (rr1_mindeltav(1), rr1_mindeltav(2), rr1_mindeltav(3), 50, 'k', 'filled')
scatter3 (rr2_mindeltav(1), rr2_mindeltav(2), rr2_mindeltav(3), 50, 'k', 'filled')
grid on
legend ('SOLE', 'Orbita Terrestre', 'Orbita Asteroide', 'Orbita di Trasferimento')

% min dt
figure_dt = figure;
scatter3 (0,0,0,100,"yellow", 'filled')
hold on
plotOrbit(O_start, 0, 2*pi, dth, 'b--');
plotOrbit(O_end, 0, 2*pi, dth, 'r--');
plotOrbit(O_best_dt, O_best_dt.th_t(1), O_best_dt.th_t(2), dth, 'g');
scatter3 (rr1_mindeltat(1), rr1_mindeltat(2), rr1_mindeltat(3), 30, 'k', 'filled')
scatter3 (rr2_mindeltat(1), rr2_mindeltat(2), rr2_mindeltat(3), 30, 'k', 'filled')
grid on
legend ('SOLE', 'Orbita Terrestre', 'Orbita Asteroide', 'Orbita di Trasferimento')
