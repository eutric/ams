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
plotOrbit(O_end,th0,thf,dth)
legend('orbita1','orbita2')
sist_can
[X,Y,Z]=sphere(1000);
k=0.5*10^4;
surf(k*X,k*Y,k*Z,"LineStyle","none","FaceColor","texturemap",CData=flip((imread(['umba.jpg']))));
xlim([-10^5;10^5])
ylim([-10^5;10^5])
zlim([-10^5;10^5])

%% prove passettini passettini
clear all
clc
close all


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
%plot orbita partenza e arrivo
th0=0;
thf=2*pi;
dth=pi/100;
plotOrbit(O_start,th0,thf,dth) % Plotto Orbita di Partenza
hold on
plotOrbit(O_end,th0,thf,dth)


[X,Y,Z]=sphere(1000);
k=0.5*10^4;
surf(k*X,k*Y,k*Z,"LineStyle","none","FaceColor","texturemap",CData=flip((imread(['earth.jpg']))));
xlim([-0.5*10^5;0.5*10^5])
ylim([-0.5*10^5;0.5*10^5])
zlim([-0.5*10^5;0.5*10^5])

%modifica piano

[delta_v_cp,om_f_cp,teta_cp,O_plan] = changeOrbitalPlane(O_start, O_end);
plotOrbit(O_plan,th0,thf,dth,'--')

%modifico anomalia pericentro
[delta_v_pp,th_i_pp,th_f_pp,th_best_pp,O_pp] = change_pericentre_arg(O_plan,O_end.om,teta_cp);
plotOrbit(O_pp,th0,thf,dth,'--')
[rr_pp,vv_pp] = par2car(O_pp, th_f_pp(2));
plot3(rr_pp(1),rr_pp(2),rr_pp(3),'*')
legend('orbita inizio','orbita fine','orbita cp', 'orbita pp')

%trasferimento bitangente/ biellittico
g={'ap','aa','pp','pa'};
min=0;
for type=g
    [delta_v1_bt, delta_v2_bt, delta_t_bt] = bitangentTransfer(O_pp, O_end, type);
    cost=delta_v1_bt+delta_v2_bt;
    if cost<min || min==0
        min=cost;
        type_min=type;
    end
end

[delta_v1_bet, delta_v2_bet, delta_v3_bet, delta_t1_bet, delta_t2_bet] = biellipticTransfer(O_pp,O_end, ra_t);









%% STRATEGIA STANDARD 
% Esempio da pericentro ad apocentro 
clear 
close all
clc

% Parametro attrattore
mu=398600;

% Orbita 1
O_start.a=24400.0;
O_start.e=0.728300;
O_start.i=0.104700;
O_start.OM=1.514000;
O_start.om=3.107000;
O_start.mu=mu;
th_start=1.665000;

% Orbita d'Arrivo
rr=[-12985.280000 3801.011400 8109.619300]';
vv=[-0.948600 -6.134000 1.356000]';
% 
[O_end,th_end] = car2par(rr,vv,mu);
% 
% [delta_v1,om_f, tetacp] = changeOrbitalPlane(O_start, O_end);
% % delta_t1 = TOF(O_start, th_start, tetacp);
% 
% [delta_v2,th_1,th_2,th_best] = change_pericentre_arg(O_start,om_f,th_start);
% % delta_t2 = TOF(O_start, tetacp, th_best);
% 
% % delta_t3 = TOF(O_start, th_best, 0);
% [delta_v11, delta_v22, delta_t4] = bitangentTransfer(O_start, O_end, 'pa');

teta1=0;
teta2=360;
th_1=teta1/180*pi;
th_2=teta2/180*pi;
delta_t5 = TOF(O_end,th_1,th_2);

% 
% 
% delta_vtot = delta_v1 + delta_v2 + delta_v11 + delta_v22;
% delta_ttot = delta_t1 + delta_t2 + delta_t3 + delta_t4 + delta_t5;
[delta_v2,th_1,th_2,th_best] = change_pericentre_arg(O_start,om_f,th_start);
delta_t2 = TOF(O_start, tetacp, th_best);

delta_t3 = TOF(O_start, th_best, 0);
[delta_v11, delta_v22, delta_t4] = bitangentTransfer(O_start, O_end, 'pa');


delta_t5 = TOF(O_end, pi/2, 3/2*pi);



delta_vtot = delta_v1 + delta_v2 + delta_v11 + delta_v22;
delta_ttot = delta_t1 + delta_t2 + delta_t3 + delta_t4 + delta_t5;

%% Appunti scenario 2
% questo giro è un trasferimento diretto
% Si ignora che i pianeti si muovono durante il traserimento (quindi la
% vera posizione dei corpi attrattori)
