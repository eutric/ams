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

% Approccio a tappeto
n = 10; % Dimensioni griglia
a = [3.4, 3.4, 1]; % Si possono cambiare ad ogni lancio per infittire dove serve
b = [3.7, 3.7, 2*pi]; 

fun_ob = @(O1, O2, x) O_tfun(O1, O2, x(1), x(2), x(3), [0,0]).cost;
[O_best, th_best_12] = tappeto(O_start, O_end, fun_ob, n, a, b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1 scenario 2, subplot di orbita terrestre 3D e piano
figure
subplot (1,2,1)
scatter3 (0,0,0, 100, [1, 0.5, 0], 'filled')
hold on
plotOrbit (O_start, 0, 2*pi, dth, 'b');
grid on
legend ('SOLE', 'Orbita terrestre')
subplot (1,2,2)
scatter (0,0,100, [1, 0.5, 0], 'filled')
hold on
sist_can(1);
plotOrbit_plane(O_start, 0, 2*pi, dth, 'blue');
legend ('SOLE', 'Orbita terrestre')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2 scenario 2 Latex - Subplot di orbita 3D e 2D

figure
subplot (1,2,1)
scatter3 (0,0,0, 100, [1, 0.5, 0], 'filled')
hold on
plotOrbit (O_end, 0, 2*pi, dth, 'r');
grid on
legend ('SOLE', 'Orbita Asteroide')
subplot (1,2,2)
scatter (0,0,100, [1, 0.5, 0], 'filled')
hold on
sist_can(1);
plotOrbit_plane(O_end, 0, 2*pi, dth, 'r');
legend ('SOLE', 'Orbita Asteroide')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
scatter3 (0,0,0, 100, [1, 0.5, 0], 'filled')
hold on
plotOrbit (O_start, 0, 2*pi, dth, 'b');
plotOrbit (O_end, 0, 2*pi, dth, 'r');
legend ('SOLE', 'Orbita Terrestre', 'Orbita Asteroide')


%% test funzione cicli annidati
clear
close all
clc
% 
[punto_1, ~] = par2car(O_start, th_best_12(1));
[punto1_t, v1_car] = par2car(O_best, th_t_best(1));
[punto2_t, v2_car] = par2car(O_best, th_t_best(2));
[punto_2, ~] = par2car(O_end, th_best_12(2));

scatter3 (punto1_t(1), punto1_t(2), punto1_t(3), 'red')
scatter3 (punto2_t(1), punto2_t(2), punto2_t(3), 'red')
scatter3 (punto_1(1), punto_1(2), punto_1(3), 'red')
scatter3 (punto_2(1), punto_2(2), punto_2(3), 'red')
plotOrbit(O_best, th_t_best(1), th_t_best(2), dth, 'g'); % Trasferimento vincente

th_best_12
costmin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 3 - trasferimento approccio a tappeto, delta v, latex
figure
scatter3 (0,0,0,100,"yellow", 'filled')
hold on
plotOrbit(O_start, 0, 2*pi, dth, 'b--');
plotOrbit(O_end, 0, 2*pi, dth, 'r--');
plotOrbit(O_best, O_best.th_t(1), O_best.th_t(2), dth, 'g');
scatter3 (punto1_t(1), punto1_t(2), punto1_t(3), 50, 'k', 'filled')
scatter3 (punto_2(1), punto_2(2), punto_2(3), 50, 'k', 'filled')
grid on
legend ('SOLE', 'Orbita Terrestre', 'Orbita Asteroide', 'Orbita di Trasferimento')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test ga
clear
close all
clc



% Terra
O_start.a = 1.4946e8;  % [km]
O_start.e = 0.016;     % [ ]
O_start.i = 9.1920e-5; % [rad]
O_start.OM = 2.7847;   % [rad]
O_start.om = 5.2643;   % [rad]
O_start.mu = mu;       % [km^3/s^2]
% r_terra = 6378.388;

[rr_start, ~] = par2car(O_start, 0);

% Asteroide 163899 (2003 SD220)
O_end.a = 0.827903 *1.496e8; % [km]
O_end.e = 0.209487;          % [ ]
O_end.i = 8.55 *pi/180;      % [rad]
O_end.OM = 273.63 *pi/180;   % [rad]
O_end.om = 327.03 *pi/180;   % [rad]
O_end.mu = mu;          

[rr_end, ~] = par2car(O_end, 0);

% fun=@(th1,th2,om)brutta(th1,th2,om,O_start,O_end);
% fun=@(x)brutta(x(1),x(2),x(3),O_start,O_end);
fun_costo = @(x)O_tfun(O_start,O_end,x(1),x(2),x(3),[0,0]).cost;
x = ga(fun_costo,3,[],[],[],[],[3.4,3.4,3.4],[3.8,3.8,3.8]);
O_opt_costo = O_tfun(O_start, O_end,x(1),x(2),x(3),[0,0])

figure
scatter3 (0, 0, 0, 'y', 'filled', LineWidth=500)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3), 'green', 'filled', LineWidth=10)
scatter3(rr_end(1), rr_end(2), rr_end(3), 'magenta', 'filled', LineWidth=10)
plotOrbit (O_start, 0, x(1), dth, 'k--');
[rr, r_theta_transfer] = plotOrbit (O_opt_costo, O_opt_costo.th_t(1), O_opt_costo.th_t(2), dth, 'r');
plotOrbit (O_end, x(2), 2*pi, dth, 'k--');
% plotOrbit (O_start, 0, 2*pi, dth, 'k--');
% plotOrbit (O_opt, 0, 2*pi, dth, 'k--');
% plotOrbit (O_end, 0, 2*pi, dth, 'k--');
% ga() trova in th1 = 3.6181, th2 = 3.4982, om = 3.7291, Dv = 6.3190, con
% durata di 11896150 s = 3304:29:10 (damn)

% O_t.a = 1.4812e+08;
% O_t.e = 0.0283;
% O_t.i = 0.0440;
% O_t.OM = 5.3829;
% O_t.om = 3.7291;
% O_t.mu = 1.327124400420000e+11;
% th_t1 = 2.5551;
% th_t2 = 4.8716;


legend ('SOLE', 'START', 'END', 'Orbita terrestre', 'Orbita di trasferimento', 'Orbita Asteroide')

O_opt_costo.tempo = seconds(O_opt_costo.tempo);
O_opt_costo.tempo.Format = 'hh:mm:ss';
O_opt_costo.tempo

%% Dati definitivi Trasferimenti
clear
close all
clc

mu_S = 1.32712440042e20 * 0.001^3; % km^3/s^2
dth = 0.001;

% Terra
O_start.a = 1.4946e8;  % [km]
O_start.e = 0.016;     % [ ]
O_start.i = 9.1920e-5; % [rad]
O_start.OM = 2.7847;   % [rad]
O_start.om = 5.2643;   % [rad]
O_start.mu = mu_S;       % [km^3/s^2]

% Asteroide 163899 (2003 SD220)
O_end.a = 0.827903 *1.496e8; % [km]
O_end.e = 0.209487;          % [ ]
O_end.i = 8.55 *pi/180;      % [rad]
O_end.OM = 273.63 *pi/180;   % [rad]
O_end.om = 327.03 *pi/180;   % [rad]
O_end.mu = mu_S;

% Min delta v
O_best_dv.a = 1.4812e+08;
O_best_dv.e = 0.0283;
O_best_dv.i = 0.0440;
O_best_dv.OM = 5.3829;
O_best_dv.om = 3.7291;
O_best_dv.mu = 1.327124400420000e+11;
O_best_dv.th_t = [2.5551, 4.8716];

% x min dv
%
% theta 1 = 3.6181  theta 2 = 3.4982  omega = 3.7291

% Min delta t
O_best_dt.a = 3.2632e+13;
O_best_dt.e = 0.9999;
O_best_dt.i = 1.5824;
O_best_dt.OM = 1.5582;
O_best_dt.om = 1.7471;
O_best_dt.mu = mu_S;
O_best_dt.th_t = [4.5360, 4.5475];

% x min dt =
% 
%     theta 1 = 6.0755   theta 2 = 3.6401   omega = 1.7471

% Punti di trasferimento
[rr1_mindeltav, vv1_mindeltav] = par2car (O_best_dv, O_best_dv.th_t(1));
[rr2_mindeltav, vv2_mindeltav] = par2car (O_best_dv, O_best_dv.th_t(2));

[rr1_mindeltat, vv1_mindeltat] = par2car (O_best_dt, O_best_dt.th_t(1));
[rr2_mindeltat, vv2_mindeltat] = par2car (O_best_dt, O_best_dt.th_t(2));
%% Scenario 2 - ottimizzo il TOF con ga()
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
r_terra = 6378.388;

[rr_start, vv_start] = par2car(O_start, 0);

% Asteroide 163899 (2003 SD220)
O_end.a = 0.827903 *1.496e8; % [km]
O_end.e = 0.209487;          % [ ]
O_end.i = 8.55 *pi/180;      % [rad]
O_end.OM = 273.63 *pi/180;   % [rad]
O_end.om = 327.03 *pi/180;   % [rad]
O_end.mu = mu;          

[rr_end, vv_end] = par2car(O_end, 0);

% fun=@(th1,th2,om)brutta(th1,th2,om,O_start,O_end);
% fun=@(x)brutta(x(1),x(2),x(3),O_start,O_end);
fun_tempo = @(x)O_tfun(O_start,O_end,x(1),x(2),x(3)).tempo;
x = ga(fun_tempo,3,[],[],[],[],[0,0,0],[2*pi,2*pi,2*pi])
O_opt_tempo = O_tfun(O_start, O_end,x(1),x(2),x(3))
O_opt_tempo.tempo = seconds(real(O_opt_tempo.tempo));
O_opt_tempo.tempo.Format = 'hh:mm:ss';
O_opt_tempo.tempo
% x =
% 
%     6.0755    3.6401    1.7471
% 
% 
% O_opt_tempo = 
% 
%   struct with fields:
% 
%         a: 3.2632e+13
%         e: 1.0000
%         i: 1.5824
%        OM: 1.5582
%        om: 1.7471
%        mu: 1.3271e+11
%     cost1: 52.2794
%     cost2: 51.0782
%      cost: 103.3576
%     tempo: 61293
%      th_t: [4.5360 4.5475]
% 
% 
% ans = 
% 
%   duration
% 
%    17:01:33


%% Scenario 2 - Ottimizzo il tempo con approccio a tappeto
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

n = 10; % Dimensioni griglia
a = 3.5;
b = 6.2;
OO_t = []; % Vettore di oggetti orbita, non si può preallocare
kk = 0; % Contatore
tmin=1e16;

figure
scatter3 (0,0,0, 'red', LineWidth=2)
hold on
plotOrbit (O_start, 0, 2*pi, dth, 'b');
plotOrbit (O_end, 0, 2*pi, dth, 'k');
tic
for th1i=linspace(5.9, 6, n)
    for th2f=linspace(3.4, 3.45, n)
        for om=linspace(0, 0.3,n)
            [O_t] = O_tfun(O_start,O_end,th1i,th2f,om);
            if O_t.e<1 && O_t.e>=0
                OO_t=[OO_t,O_t];
                if O_t.tempo<tmin
                    tmin=O_t.tempo;
                    O_best=O_t;
                    th_t_best=O_t.th_t;
                    th_best_12 = [th1i, th2f];
                end
            else 
                kk = kk+1;
            end

        end
    end
end
t_elaborazione = toc 
[punto_1, ~] = par2car(O_start, th_best_12(1));
[punto1_t, ~] = par2car(O_best, th_t_best(1));
[punto2_t, ~] = par2car(O_best, th_t_best(2));
[punto_2, ~] = par2car(O_end, th_best_12(2));

scatter3 (punto1_t(1), punto1_t(2), punto1_t(3), 'red')
scatter3 (punto2_t(1), punto2_t(2), punto2_t(3), 'red')
scatter3 (punto_1(1), punto_1(2), punto_1(3), 'red')
scatter3 (punto_2(1), punto_2(2), punto_2(3), 'red')
plotOrbit(O_best, th_t_best(1), th_t_best(2), dth, 'g'); % Trasferimento vincente

O_best.tempo = seconds(O_best.tempo);
O_best.tempo.Format = 'hh:mm:ss';
O_best.tempo
th_best_12
O_best.om
file_O_best = fopen ('O_best.txt', 'a');
fprintf (file_O_best, 'Risultato - %s\n', string(datetime('now')));
fprintf (file_O_best, 'In %d secondi, con n = %d\n', t_elaborazione, n);
fprintf (file_O_best, 'La miglior orbita in termini di tempo, trovata con cicli annidati ha\n');
fprintf (file_O_best, 'a = %f\n', O_best.a);
fprintf (file_O_best, 'e = %f\n', O_best.e);
fprintf (file_O_best, 'i = %f\n', O_best.i);
fprintf (file_O_best, 'OM = %f\n', O_best.OM);
fprintf (file_O_best, 'om = %f\n', O_best.om);
fprintf (file_O_best, 'In termini di tempo il costo è %s\n', string(O_best.tempo));
fprintf (file_O_best, 'In termini di velocità il costo è %f\n', O_best.cost);
fprintf (file_O_best, 'Il trasferimento avviene in theta_1_t = %d, theta_2_t = %d\n', O_best.th_t(1), O_best.th_t(2));
fprintf (file_O_best, 'In riferimento a alla orbita 1 e 2, %d e %d\n\n', th_best_12(1), th_best_12(2));
fclose (file_O_best);
% Con n = 50, a = 3.5, b = 6.2, per om 0-b, si trova in th1=6 .0898,
% th2=3.6653 e om=1.8980 21:03:12, in 470 secondi di computazione

%% 
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


O_best.a = 974121862.949664;
O_best.e = 0.945997;
O_best.i = 0.633663;
O_best.OM = 1.572542;
O_best.om = 1.897959;
O_best.mu = mu;
th_t_best = [4.385081, 4.398208];
th_best_12 = [6.089796, 3.665306];
[punto_1, ~] = par2car(O_start, th_best_12(1));
[punto1_t, ~] = par2car(O_best, th_t_best(1));
[punto2_t, ~] = par2car(O_best, th_t_best(2));
[punto_2, ~] = par2car(O_end, th_best_12(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 4 - scenario 2, trasferimento approccio a tappeto min (delta t)
figure
scatter3 (0,0,0,100,"yellow", 'filled')
hold on
plotOrbit(O_start, 0, 2*pi, dth, 'b--');
plotOrbit(O_end, 0, 2*pi, dth, 'r--');
plotOrbit(O_best, th_t_best(1), th_t_best(2), dth, 'g');
scatter3 (punto1_t(1), punto1_t(2), punto1_t(3), 30, 'k', 'filled')
scatter3 (punto_2(1), punto_2(2), punto_2(3), 30, 'k', 'filled')
grid on
legend ('SOLE', 'Orbita Terrestre', 'Orbita Asteroide', 'Orbita di Trasferimento')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ottimizzazione pesata
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
r_terra = 6378.388;

[rr_start, vv_start] = par2car(O_start, 0);

% Asteroide 163899 (2003 SD220)
O_end.a = 0.827903 *1.496e8; % [km]
O_end.e = 0.209487;          % [ ]
O_end.i = 8.55 *pi/180;      % [rad]
O_end.OM = 273.63 *pi/180;   % [rad]
O_end.om = 327.03 *pi/180;   % [rad]
O_end.mu = mu;          

[rr_end, vv_end] = par2car(O_end, 0);

% fun=@(th1,th2,om)brutta(th1,th2,om,O_start,O_end);
% fun=@(x)brutta(x(1),x(2),x(3),O_start,O_end);
fun_opt = @(x)O_tfun(O_start,O_end,x(1),x(2),x(3),[50,1]).fun_pesata;
x = ga(fun_opt,3,[],[],[],[],[0,0,0],[2*pi,2*pi,2*pi])
O_opt_pesata = O_tfun(O_start, O_end,x(1),x(2),x(3), [0,0])
O_opt_pesata.tempo = seconds(O_opt_pesata.tempo);
O_opt_pesata.tempo.Format = 'hh:mm:ss';
O_opt_pesata.tempo
O_opt_pesata.cost

% alpha = 50, beta = 1
% x =
% 
%     6.0280    3.6867    2.1148
% 
% 
% O_opt_pesata = 
% 
%   struct with fields:
% 
%              a: 1.3704e+08
%              e: 0.2160
%              i: 0.0491
%             OM: 1.5125
%             om: 2.1148
%             mu: 1.3271e+11
%          cost1: 6.2575
%          cost2: 6.2823
%           cost: 12.5399
%          tempo: 4.7740e+05
%     fun_pesata: 0
%           th_t: [4.1666 4.2602]
% 
% 
% ans = 
% 
%   duration
% 
%    132:36:35
% 
% 
% ans =
% 
%    12.5399

O_best.a = 1.3704e+08;
O_best.e = 0.2160;
O_best.i = 0.0491;
O_best.OM = 1.5125;
O_best.om = 2.1148;
O_best.cost = 12.5399;
O_best.tempo = 4.7740e05;
O_best.th_t = [4.1666, 4.2602];
O_best.mu = mu;
th_best_12 = [6.0280, 3.6867];

[punto1_t, ~] = par2car(O_best, O_best.th_t(1));
[punto_2, ~] = par2car(O_end, th_best_12(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 4 - scenario 2, trasferimento approccio a tappeto min (delta t)
figure
scatter3 (0,0,0,100,"yellow", 'filled')
hold on
plotOrbit(O_start, 0, 2*pi, dth, 'b--');
plotOrbit(O_end, 0, 2*pi, dth, 'r--');
plotOrbit(O_best, O_best.th_t(1), O_best.th_t(2), dth, 'g');
scatter3 (punto1_t(1), punto1_t(2), punto1_t(3), 30, 'k', 'filled')
scatter3 (punto_2(1), punto_2(2), punto_2(3), 30, 'k', 'filled')
grid on
legend ('SOLE', 'Orbita Terrestre', 'Orbita Asteroide', 'Orbita di Trasferimento')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% 
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


O_best.a = 1.3535e+08;
O_best.e = 0.1959;
O_best.i = 0.0272;
O_best.OM = 1.5022;
O_best.om = 2.1982;
O_best.mu = mu;
th_t_best = [4.0817, 4.1967];
th_best_12 = [6.0164, 3.6966];
[punto_1, ~] = par2car(O_start, th_best_12(1));
[punto1_t, ~] = par2car(O_best, th_t_best(1));
[punto2_t, ~] = par2car(O_best, th_t_best(2));
[punto_2, ~] = par2car(O_end, th_best_12(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 4 - scenario 2, trasferimento approccio a tappeto min (delta t)
figure
scatter3 (0,0,0,100,"yellow", 'filled')
hold on
plotOrbit(O_start, 0, 2*pi, dth, 'b--');
plotOrbit(O_end, 0, 2*pi, dth, 'r--');
plotOrbit(O_best, th_t_best(1), th_t_best(2), dth, 'g');
scatter3 (punto1_t(1), punto1_t(2), punto1_t(3), 30, 'k', 'filled')
scatter3 (punto_2(1), punto_2(2), punto_2(3), 30, 'k', 'filled')
grid on
legend ('SOLE', 'Orbita Terrestre', 'Orbita Asteroide', 'Orbita di Trasferimento')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fmincon per il deltav
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
% r_terra = 6378.388;

[rr_start, ~] = par2car(O_start, 0);

% Asteroide 163899 (2003 SD220)
O_end.a = 0.827903 *1.496e8; % [km]
O_end.e = 0.209487;          % [ ]
O_end.i = 8.55 *pi/180;      % [rad]
O_end.OM = 273.63 *pi/180;   % [rad]
O_end.om = 327.03 *pi/180;   % [rad]
O_end.mu = mu;          

[rr_end, ~] = par2car(O_end, 0);

% fun=@(th1,th2,om)brutta(th1,th2,om,O_start,O_end);
% fun=@(x)brutta(x(1),x(2),x(3),O_start,O_end);
fun_costo = @(x)O_tfun(O_start,O_end,x(1),x(2),x(3),[0,0]).cost;
x0 = [pi,pi,pi];
x = fmincon(fun_costo,x0,[],[],[],[],[0,0,0],[2*pi, 2*pi, 2*pi])
O_opt_costo = O_tfun(O_start, O_end,x(1),x(2),x(3),[0,0])

figure
scatter3 (0, 0, 0, 'y', 'filled', LineWidth=500)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3), 'green', 'filled', LineWidth=10)
scatter3(rr_end(1), rr_end(2), rr_end(3), 'magenta', 'filled', LineWidth=10)
plotOrbit (O_start, 0, x(1), dth, 'k--');
[rr, r_theta_transfer] = plotOrbit (O_opt_costo, O_opt_costo.th_t(1), O_opt_costo.th_t(2), dth, 'r');
plotOrbit (O_end, x(2), 2*pi, dth, 'k--');
% plotOrbit (O_start, 0, 2*pi, dth, 'k--');
% plotOrbit (O_opt, 0, 2*pi, dth, 'k--');
% plotOrbit (O_end, 0, 2*pi, dth, 'k--');
% ga() trova in th1 = 3.6181, th2 = 3.4982, om = 3.7291, Dv = 6.3190, con
% durata di 3304:29:10 (damn)

% O_t.a = 1.4812e+08;
% O_t.e = 0.0283;
% O_t.i = 0.0440;
% O_t.OM = 5.3829;
% O_t.om = 3.7291;
% O_t.mu = 1.327124400420000e+11;
% th_t1 = 2.5551;
% th_t2 = 4.8716;


legend ('SOLE', 'START', 'END', 'Orbita terrestre', 'Orbita di trasferimento', 'Orbita Asteroide')

O_opt_costo.tempo = seconds(O_opt_costo.tempo);
O_opt_costo.tempo.Format = 'hh:mm:ss';
O_opt_costo.tempo

%% fmincon per il deltat
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
% r_terra = 6378.388;

[rr_start, ~] = par2car(O_start, 0);

% Asteroide 163899 (2003 SD220)
O_end.a = 0.827903 *1.496e8; % [km]
O_end.e = 0.209487;          % [ ]
O_end.i = 8.55 *pi/180;      % [rad]
O_end.OM = 273.63 *pi/180;   % [rad]
O_end.om = 327.03 *pi/180;   % [rad]
O_end.mu = mu;          

[rr_end, ~] = par2car(O_end, 0);

% fun=@(th1,th2,om)brutta(th1,th2,om,O_start,O_end);
% fun=@(x)brutta(x(1),x(2),x(3),O_start,O_end);
fun_tempo = @(x)O_tfun(O_start,O_end,x(1),x(2),x(3),[0,0]).tempo;
x0 = [6, 3.6, 1.7]
x = fmincon(fun_tempo,x0,[],[],[],[],[0,0,0],[2*pi, 2*pi, 2*pi])
O_opt_costo = O_tfun(O_start, O_end,x(1),x(2),x(3),[0,0])

figure
scatter3 (0, 0, 0, 'y', 'filled', LineWidth=500)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3), 'green', 'filled', LineWidth=10)
scatter3(rr_end(1), rr_end(2), rr_end(3), 'magenta', 'filled', LineWidth=10)
plotOrbit (O_start, 0, x(1), dth, 'k--');
[rr, r_theta_transfer] = plotOrbit (O_opt_costo, O_opt_costo.th_t(1), O_opt_costo.th_t(2), dth, 'r');
plotOrbit (O_end, x(2), 2*pi, dth, 'k--');
% plotOrbit (O_start, 0, 2*pi, dth, 'k--');
% plotOrbit (O_opt, 0, 2*pi, dth, 'k--');
% plotOrbit (O_end, 0, 2*pi, dth, 'k--');
% ga() trova in th1 = 3.6181, th2 = 3.4982, om = 3.7291, Dv = 6.3190, con
% durata di 3304:29:10 (damn)

% O_t.a = 1.4812e+08;
% O_t.e = 0.0283;
% O_t.i = 0.0440;
% O_t.OM = 5.3829;
% O_t.om = 3.7291;
% O_t.mu = 1.327124400420000e+11;
% th_t1 = 2.5551;
% th_t2 = 4.8716;


legend ('SOLE', 'START', 'END', 'Orbita terrestre', 'Orbita di trasferimento', 'Orbita Asteroide')

O_opt_costo.tempo = seconds(O_opt_costo.tempo);
O_opt_costo.tempo.Format = 'hh:mm:ss';
O_opt_costo.tempo