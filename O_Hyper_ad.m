function [O_hyper_Ad,DELTA_V_AD] = O_Hyper_ad(theta_parc)
% Dati

G = 6.67e-11; % N*m^2/kg^2

m_T = 5.974e24;  % kg
m_S = 1.989e30;  % kg
m_A = 5.1827e11; % Kg

R_T = 6378;   % km
R_S = 696340; % km
R_A = 0.79/2; % km

mu_T = 398600;  % km^3/s^2
mu_S = 1.32712440042e20 * 0.001^3; % km^3/s^2
mu_A = G * m_A;


dth = 0.1;

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

[rr_t1, vv_t1, matrici] = par2car(O_t_scenario2, th1_transfer);
vv_inf1 = vv_t1 - vv_end1_S;
%Prova iperbole con modifica di piano
%1)modifico sistema riferimento passando da eclittico(sole) a eci(terra)
T_eci_eclip=[1 0 0;
            0 cos(eps) sin(eps);
            0 -sin(eps) cos(eps)];
vv_inf1_eci=T_eci_eclip'*vv_inf1;
rr_inf_eci=vv_inf1_eci/norm(vv_inf1_eci);

COST_vect=[];
cost_min=10;
%2)scelgo r_h per esempio raggio di 90 gradi su orbita di parcheggio dell'orbita di parcheggio
[rr_h_earth,vv_h_earth]=par2car(O_parcheggio,theta_parc);
rr_h=rr_h_earth;
r_h=norm(rr_h);
alpha=acos(dot(rr_h,rr_inf_eci)./norm(rr_h)./norm(rr_inf_eci));
a_h=-O_parcheggio.mu/(norm(vv_inf1_eci).^2);
%X=[eh,theta_inf,theta_h]
fun=@(X)[(a_h*(1-X(1).^2))./(1+X(1).*cos(X(3)))-r_h;
    cos(X(2))+1./X(1);
    X(3)+alpha-X(2)
    ];
x0=[1,pi/2,0];
X=fsolve(fun,x0);
eh=X(1);
theta_inf=X(2);
theta_h=X(3);
if theta_inf>pi/2 && theta_inf>pi
    printf('modificare il raggio di intersezione')
    stop
end

%definisco il piano dell'iperbole uguale al piano dell orbita di
%trasferimento delo sc2
hh_hyper=cross(rr_h,rr_inf_eci)/norm(cross(rr_inf_eci,rr_h));
NN_h=cross([0;0;1],hh_hyper)/norm(cross([0;0;1],hh_hyper));
O_hyper_Ad.i=acos(dot(hh_hyper,[0,0,1]));
if NN_h(2)>=0
    O_hyper_Ad.OM=acos(NN_h(1));
elseif NN_h(2)<0
     O_hyper_Ad.OM=2*pi-acos(NN_h(1));
end
%uso sistema
fun2=@(e)[dot(rr_h,e)./norm(rr_h)-cos(theta_h);
    dot(rr_inf_eci,e)-cos(theta_inf);
    norm(e)-1
    ];
e0=[1,1,1];
ee_h=fsolve(fun2,e0);


if dot(ee_h,[0,0,1])>=0
    O_hyper_Ad.om=acos(dot(NN_h,ee_h));
elseif dot(ee_h,[0,0,1])<0
    O_hyper_Ad.om=2*pi-acos(dot(NN_h,ee_h));
end


O_hyper_Ad.a=-mu_T/norm(vv_inf1).^2;
O_hyper_Ad.e=eh;


%calcolo costi%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O_hyper_Ad.mu=mu_T;
[rr_hc,vv_h]=par2car(O_hyper_Ad,-theta_h);%uso - angolo perchè il vettore rr_h si trova prima del pericentro
[rr_pc,vv_p]=par2car(O_parcheggio,theta_parc);
DELTA_V_AD=norm(vv_h-vv_p);
end

