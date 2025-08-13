clear
close all
clc

m_T = 5.974e24; % kg
m_S = 1.989e30; % kg

mu_T = 398600;  % km^3/s^2
mu_S = 1.32712440042e20 * 0.001^3; % km^3/s^2

R_T = 6378; % km
R_S = 696340; % km
R_t_s = 149597870707 * 10^-3; % km

% !!!! RIFERITA ALLA TERRA
O_parcheggio.a = 3.918174284575137e+04; % km
O_parcheggio.e = 0.597404804289614;
O_parcheggio.i = 0.590342952766537;
O_parcheggio.OM = 1.749480463333461;
O_parcheggio.om = 1.177309652253865;
O_parcheggio.mu = 398600;
[rr_start_T, vv_start_T] = par2car(O_parcheggio, 0); % Per ora mi metto nel pericentro; SdR Terra

% Sarebbe l'orbita di parcheggio del nostro gruppo, dello scenario 1

% Terra !!!! RIFERITA AL SOLE
O_start_scenario2.a = 1.4946e8;  % [km]
O_start_scenario2.e = 0.016;     % [ ]
O_start_scenario2.i = 9.1920e-5; % [rad]
O_start_scenario2.OM = 2.7847;   % [rad]
O_start_scenario2.om = 5.2643;   % [rad]
O_start_scenario2.mu = mu_S;       % [km^3/s^2]
th1_scenario2 = 3.618362149407213; % [rad]
[rr_end1_S, vv_end1_S] = par2car (O_start_scenario2, th1_scenario2); % SdR Sole

% Trasferimento dello scenario 2 !!! Riferita al sole
O_t_scenario2.a = 1.481231188816134e+08;
O_t_scenario2.e = 0.028295312693238;
O_t_scenario2.i = 0.043993770239157;
O_t_scenario2.OM = 5.382908205403476;
O_t_scenario2.om = 3.729157082348018;
O_t_scenario2.mu = 1.327124400420000e+11;
th1_transfer = 2.555106932985983;
th2_transfer = 4.871467598525006;

[rr_t1, vv_t1] = par2car(O_t_scenario2, th1_transfer); % Coordinate cartesiane 
% rispetto al sole del satellite quando fa la manovra, la velocità è la
% velocità subito dopo l'impulso, quindi sta già percorrendo l'orbita di
% trasferimento

% Per approssimazione della patched conics, uscito dalla SOI terrestre,
% posso confondere la posizione del satellite con quella del pianeta,
% quindi, centrandomi nella terra, scelgo il punto in cui passo ad
% un'iperbole e da lì vado all'infinito
% Per questa prima iperbole il piano è lo stesso dell'orbita di parcheggio,
% visto che bisogna aumentare l'energia, direi che conviene fissare il
% pericentro con il pericentro dell'orbita di parcheggio
% Voglio arrivare all'infinito con, nel SdR Sole, la velocità che mi
% permette di percorre l'orbita di trasferimento dello scenario 2, è così
% imposta la v_inf dell'iperbole rispetto alla terra

vv_inf1 = vv_t1 - vv_end1_S; % Velocità dell'orbita di T - la velocità 
% della terra ==> Velocità infinito dell'iperbole relativa alla terra
% ! Dovrebbe essere rispetto alla terra.

% Il raggio di pericentro dell'iperbole, lo facciamo coincidere col
% pericentro del'orbita di parcheggio:
rp_parcheggio = r_parametrica(O_parcheggio, 0); % rp_iperbole1
vp_parcheggio = v_theta(O_parcheggio, 0);

% Caratterizzo l'orbita iperbolica, per trovarne v nel pericentro
O_hyper1.a = -mu_T/norm(vv_inf1)^2; % Dalla conservazione dell'energia meccanica

% Passando per l'equazione della conica al pericentro, si trova un'eq di
% secondo grado per le eccentricità
eH_1 = (-rp_parcheggio/2 - ...
    sqrt( (rp_parcheggio/2)^2 - O_hyper1.a*(rp_parcheggio-O_hyper1.a) ))/O_hyper1.a; % 1.0858 > 1 ==> ok
eH_2 = (-rp_parcheggio/2 + ...
    sqrt( (rp_parcheggio/2)^2 - O_hyper1.a*(rp_parcheggio-O_hyper1.a) ))/O_hyper1.a; % < 0 la scarto

O_hyper1.e = eH_1;
O_hyper1.i = 0.590342952766537;
O_hyper1.OM = 1.749480463333461;
O_hyper1.om = 1.177309652253865; % Stesso piano e orientazione dell'orbita di parcheggio
O_hyper1.mu = 398600;

% Trovo la velocità nel pericentro dell'orbita:
vp_hyper = v_theta(O_hyper1, 0)

% e quindi il deltaV richiesto, visto che il sistema di riferimento è lo
% stesso, ho preso le parametriche, toccherebbe verificare che sia giusto,
% ma devo andare
DELTA_V1 = abs(vp_hyper - vp_parcheggio) % Bello basso
%%%%%%%%%%%%%%%%%%%%%%% Parentesi Grafica %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R_SOI_TERRA = R_t_s* (m_T/m_S)^(2/5) % km
% 
% [X, Y, Z] = sphere;
% X_SOI_TERRA = R_SOI_TERRA * X;
% Y_SOI_TERRA = R_SOI_TERRA * Y;
% Z_SOI_TERRA = R_SOI_TERRA * Z;
% 
% X_T = R_T * X;
% Y_T = R_T * Y;
% Z_T = R_T * Z;
% 
% figure
% 
% surf_T = surf (X_T, Y_T, Z_T);
% shading interp
% surf_T.FaceAlpha = 0.3;
% surf_T.EdgeColor = 'none';
% 
% hold on
% 
% surf_SOI_TERRA = surf (X_SOI_TERRA, Y_SOI_TERRA, Z_SOI_TERRA);
% shading interp
% surf_SOI_TERRA.FaceAlpha = 0.3;
% surf_SOI_TERRA.EdgeColor = 'none';
% 
% scatter3(rr_start(1), rr_start(2), rr_start(3), 'filled', LineWidth=1)
% scatter3(rr_end(1), rr_end(2), rr_end(3), 'filled', LineWidth=1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%