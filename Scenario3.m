clear
close all
clc

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
[rr_start_T, vv_start_T] = par2car(O_parcheggio, 6.2832); % Per ora mi metto nel pericentro; SdR Terra


% Sarebbe l'orbita di parcheggio del nostro gruppo, dello scenario 1

% Terra RIFERITA AL SOLE
O_start_scenario2.a = 1.4946e8;  % [km]
O_start_scenario2.e = 0.016;     % [ ]
O_start_scenario2.i = 9.1920e-5; % [rad]
O_start_scenario2.OM = 2.7847;   % [rad]
O_start_scenario2.om = 5.2643;   % [rad]
O_start_scenario2.mu = mu_S;       % [km^3/s^2]
th1_scenario2 = 3.618362149407213; % [rad]
[rr_end1_S, vv_end1_S] = par2car (O_start_scenario2, th1_scenario2); % SdR Sole

% Trasferimento dello scenario 2  Riferita al sole
O_t_scenario2.a = 1.481233898202188e+08;
O_t_scenario2.e = 0.028294269485414;
O_t_scenario2.i = 0.043999052217133;
O_t_scenario2.OM = 5.382840011661747;
O_t_scenario2.om = 3.729187666090036;
O_t_scenario2.mu = 1.327124400420000e+11;
th1_transfer = 2.555076341906757;
th2_transfer = 4.871489330787082;

[rr_t1, vv_t1, matrici] = par2car(O_t_scenario2, th1_transfer); % Coordinate cartesiane 
% rispetto al sole del satellite quando fa la manovra, la velocità è la
% velocità subito dopo l'impulso, quindi sta già percorrendo l'orbita di
% trasferimento


vv_inf1 = vv_t1 - vv_end1_S; % Velocità dell'orbita di T - la velocità 
% della terra ==> Velocità infinito dell'iperbole relativa alla terra
% Dovrebbe essere rispetto alla terra.

% Il raggio di pericentro dell'iperbole, lo facciamo coincidere col
% pericentro del'orbita di parcheggio:
rp_parcheggio = r_parametrica(O_parcheggio, 0); % rp_iperbole1
vp_parcheggio = v_theta(O_parcheggio, 0)

% Caratterizzo l'orbita iperbolica, per trovarne v nel pericentro
O_hyper1.a = -mu_T/norm(vv_inf1)^2; % Dalla conservazione dell'energia meccanica

% Passando per l'equazione della conica al pericentro
O_hyper1.e = (O_hyper1.a-rp_parcheggio)/O_hyper1.a;
O_hyper1.i = 0.590342952766537;
O_hyper1.OM = 1.749480463333461;
O_hyper1.om = 1.177309652253865; % Stesso piano e orientazione dell'orbita di parcheggio
O_hyper1.thetainf = acos(-1/O_hyper1.e); 
O_hyper1.mu = 398600;
O_hyper1.delta=2*acos(1/O_hyper1.e);
O_hyper1.Delta=-O_hyper1.a*O_hyper1.e*cos(O_hyper1.delta/2);

% Trovo la velocità nel pericentro dell'orbita:
vp_hyper = v_theta(O_hyper1, 0)

% e quindi il deltaV richiesto, visto che il sistema di riferimento è lo
% stesso, ho preso le parametriche, toccherebbe verificare che sia giusto,
% ma devo andare
DELTA_V1 = abs(vp_hyper - vp_parcheggio)

%COSTO in tempo: prendo la legge oraria dell iperbole artendo da un th = 0
%(pericentro iperbole, arrivando a th tale da raggiungere il raggio della
%sfera di influenza
d_ear_sun=norm(rr_end1_S);
r_soi_ear=d_ear_sun*(m_T/m_S)^(2/5);
th_r_soi_ear=acos(1/O_hyper1.e*(O_hyper1.a*(1-O_hyper1.e^2)/r_soi_ear-1));
DELTA_T_SOI1 = TOF_open(O_hyper1, 0, th_r_soi_ear)


% Vediamo se funziona, disegnando
Terra_3D(R_T)
hold on
plotOrbit (O_parcheggio, 0, 2*pi, dth, 'k--');
plotOrbit (O_hyper1, 0, O_hyper1.thetainf, dth, 'r');
xlim([-1e5,1e5]);
ylim([-1e5,1e5]);
zlim([-1e5,1e5]);

% MANOVRA DI RIENTRO IN ASTEROIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%definisco un orbita intorno all asteroide tale da essere contenuta nella
%sfera di influenza dell asteroide rispetto al sole e che sia sullo stesso
%piano di quella di traferimento nello scenario 2

% Orbita asteroide intorno al sole
O_A.a = 0.827903 *1.496e8; % [km]
O_A.e = 0.209487;          % [ ]
O_A.i = 8.55 *pi/180;      % [rad]
O_A.OM = 273.63 *pi/180;   % [rad]
O_A.om = 327.03 *pi/180;
O_A.mu=1.32712440042e20 * 0.001^3;
% th_sc2_ast= 3.5105; % Dove era stato preso?
th_sc2_ast = 3.498164148433855; % th2 trovato con ga()

%
% calcolo sfera di influenza dell asteroide rispetto al sole
% prima cosa calcolo distanza asteroide sole nel punto in cui avviene la
% manovra
d_ast_sun=r_parametrica(O_A, th_sc2_ast); % Distanza asteroide - sole al 
% momento d'arrivo
R_SOI_A=d_ast_sun*(m_A/m_S)^(2/5);



% La velocità all'infinito dell'iperbole di rientro è definita, è quella
% dell'orbita di trasferimento in theta 2
[rr_th2_t, vv_th2_t] = par2car (O_t_scenario2, th2_transfer); % Velocità rispetto al Sole
% del satellite nell'orbita di trasferimento

[rr_th2_S, vv_th2_S] = par2car (O_A, th_sc2_ast); % Velocità rispetto al Sole
% dell'asteroide, nel punto di trasferimento

vv_inf2_S = -vv_th2_t + vv_th2_S; % km/s, velocità di eccesso iperbolico, rispetto
% all'asteroide... su che piano sono?
v_inf2 = norm(vv_inf2_S);
rr_inf2_S = rr_th2_t - rr_th2_S; % km, distanza relativa rispetto allo 
% asteroide, per approssimazione delle patched conics, = 0
O_hyper2.i = O_t_scenario2.i;
O_hyper2.om = O_t_scenario2.om;
O_hyper2.OM = O_t_scenario2.OM;

O_hyper2.a = - mu_A/v_inf2^2;
O_hyper2.mu = mu_A;
% Ignorando per ora i piani in cui giaccio, posso trovare il costo minimo
% della manovra, costruendo l'ellisse con la velocità più alta possibile al
% pericentro

% bisogna quindi scegliere un orbita con raggio massimo di 5.4226 km 
%OSS: datemi comferma ma volendo ottimizzare la manovra forse conviene un
%orbita molto ellittica così la velocità di pericentro è alta e risparmio
%carburante, giac d'accordo
r_a = 0.9*R_SOI_A; % Margine del 10%, aumentandolo, aumentano i costi della
% manovra
r_p = 1.2*R_A;

O_hyper2.e = (O_hyper2.a-r_p)/O_hyper2.a;
O_hyper2.thetainf = acos(-1/O_hyper2.e);
O_hyper2.delta=2*acos(1/O_hyper2.e);
O_hyper2.Delta=-O_hyper2.a*O_hyper2.e*cos(O_hyper2.delta/2);

% Dell'orbita d'arrivo impongo solo r_a
O_arrivo.mu = mu_A;
O_arrivo.a=(r_a+r_p)/2;
O_arrivo.e=(r_a-r_p)/2/O_arrivo.a;
O_arrivo.i = O_hyper2.i;
O_arrivo.om = O_hyper2.om;
O_arrivo.OM = O_hyper2.OM;
% Quindi, l'impulso finale per pormi nell'ultima orbita di parcheggio, sarà
DELTA_V2 = abs (v_theta (O_arrivo, 0) - v_theta (O_hyper2, 0))


%COSTO in tempo: prendo la legge oraria dell iperbole artendo da un th = 0
%(pericentro iperbole, arrivando a th tale da raggiungere il raggio della
%sfera di influenza

th_r_soi_ast=acos(1/O_hyper2.e*(O_hyper2.a*(1-O_hyper2.e^2)/R_SOI_A-1));
DELTA_T_SOI2= TOF_open(O_hyper2, 0, th_r_soi_ast)

Terra_3D(0.36)
hold on
plotOrbit (O_arrivo, 0, 2*pi, dth, 'k--');
plotOrbit (O_hyper2, 0, O_hyper2.thetainf, dth, 'r');
xlim([-10,10]);
ylim([-10,10]);
zlim([-10,10]);



% AVANZATO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Prova iperbole con modifica di piano
%1)modifico sistema riferimento passando da eclittico(sole) a eci(terra)
eps=23.45*pi/180;
T_eci_eclip=[1 0 0;
            0 cos(eps) sin(eps);
            0 -sin(eps) cos(eps)];
vv_inf1_eci=T_eci_eclip'*vv_inf1;
rr_inf_eci=vv_inf1_eci/norm(vv_inf1_eci);



%2)scelgo r_h per esempio raggio di 90 gradi su orbita di parcheggio dell'orbita di parcheggio

%uso una funzione che itera i passaggi successivi per trovare quale th è il
%migliore, per poi usare questo th per il calcolo esplicito dell'iperbole
costmin=10;
for th=linspace(0,2*pi,100)
    [O_hyper_Ad,cost] = O_Hyper_ad(th);
    if cost<costmin
        costmin=cost;
        O_opt=O_hyper_Ad;
        th_opt=th;
    end
end

theta_parc=th_opt;
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
O_hyper_Ad.delta=2*acos(1/O_hyper_Ad.e);
O_hyper_Ad.Delta=-O_hyper_Ad.a*O_hyper_Ad.e*cos(O_hyper_Ad.delta/2);


%calcolo costi velocità%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O_hyper_Ad.mu=mu_T;
[rr_hc,vv_h]=par2car(O_hyper_Ad,-theta_h);%uso - angolo perchè il vettore rr_h si trova prima del pericentro
[rr_pc,vv_p]=par2car(O_parcheggio,theta_parc);
DELTA_V_AD=norm(vv_h-vv_p);




%calcolo costi tempo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_ear_sun=norm(rr_end1_S);
r_soi_ear=d_ear_sun*(m_T/m_S)^(2/5);
th_r_soi_ear2=acos(1/O_hyper_Ad.e*(O_hyper_Ad.a*(1-O_hyper_Ad.e^2)/r_soi_ear-1));
t1=TOF_open(O_hyper_Ad,theta_h,th_r_soi_ear2);

%calcolo inoltre il tempo necessario a passare da th0=0 fino a theta
%optimal

t_parch=TOF(O_parcheggio,0,th_opt);

t_tot=t_parch+t1
%plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Terra_3D(R_T)
hold on
plotOrbit(O_hyper_Ad,-theta_inf,theta_h,0.01,'r -- ');
plotOrbit(O_hyper_Ad,theta_h,theta_inf,0.01,'r');
plotOrbit(O_parcheggio,0,th_opt,0.01,'k');
plotOrbit(O_parcheggio,th_opt,2*pi,0.01,'k--');
quiver3(0,0,0,rr_h(1),rr_h(2),rr_h(3))

ee_h_plot=0.5e5*ee_h;
quiver3(0,0,0,ee_h_plot(1),ee_h_plot(2),ee_h_plot(3))
quiver3(0,0,0,1e5*hh_hyper(1),1e5*hh_hyper(2),1e5*hh_hyper(3))
quiver3(0,0,0,1e5*NN_h(1),1e5*NN_h(2),1e5*NN_h(3))
quiver3(0,0,0,1e5*rr_inf_eci(1),1e5*rr_inf_eci(2),1e5*rr_inf_eci(3))
legend('attrattore','','iperbole','orbita di partenza','','raggio di intersezione','direz. eccentricità','direz. vettore momento della qt. moto','direz. asse dei nodi','direz. velocità asintotica')
xlim([-1e5,1e5]);
ylim([-1e5,1e5]);
zlim([-1e5,1e5]);




