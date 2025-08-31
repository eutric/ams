%% Scenario 1
% Trasferimento geocentrico - da GTO ad orbita di parcheggio assegnata
clear all
close all
clc

% parametro gravitazionale della Terra
mu=398600;          % [km^3/s^2]
% raggio Terra 
r_terra = 6378.388; % [km]
% passo
dth=0.001;          % [s]

% Orbita di partenza
O_start.a = 24400.0;     % [km]
O_start.e = 0.728300;    % [ ]
O_start.i = 0.104700;    % [rad]
O_start.OM = 1.514000;   % [rad]
O_start.om = 3.107000;   % [rad]
O_start.mu = mu;         % [km^3/s^2]
th_start = 1.665000;     % [rad]

[rr_start, vv_start] = par2car(O_start, th_start);

% Orbita d'arrivo
rr_end=[-12985.280000 3801.011400 8109.619300]';
vv_end=[-0.948600 -6.134000 1.356000]';

[O_end,th_end] = car2par(rr_end,vv_end,mu);

%% Altri dati rilevanti 

% Orbita iniziale 
rp1 = O_start.a*(1-O_start.e);  % raggio di pericentro                 [km]
ra1 = O_start.a*(1+O_start.e);  % raggio di apocentro                  [km]
p1 = O_start.a*(1-O_start.e^2); % semilato retto                       [km]
E1 = 0.5 * norm(vv_start)^2 - mu/norm(rr_start); % energia specifica   [J/kg]
T1 = 2*pi * sqrt(O_start.a^3/mu); % periodo orbitale                   [s]
b1 = O_start.a*sqrt(1-O_start.e^2); % semiasse minore                  [km]

% Orbita iniziale 
rp2 = O_end.a*(1-O_end.e);  % raggio di pericentro                 [km]
ra2 = O_end.a*(1+O_end.e);  % raggio di apocentro                  [km]
p2 = O_end.a*(1-O_end.e^2); % semilato retto                       [km]
E2 = 0.5 * norm(vv_end)^2 - mu/norm(rr_end); % energia specifica   [J/kg]
T2 = 2*pi * sqrt(O_end.a^3/mu); % periodo orbitale                 [s]
b2 = O_end.a*sqrt(1-O_end.e^2); % semiasse minore                  [km]

%% Rappresentazione grafica orbita iniziale in ECI
Terra_3D(r_terra)
set(gcf, 'Name', 'Orbita di partenza', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
plotOrbit(O_start,0,2*pi,dth,'b');

%% Rappresentazione grafica orbita iniziale nel sistema perifocale
x1 = rp1;
y1 = 0;
x2 = -ra1;
y2 = 0;
t = linspace(0,2*pi);
X = O_start.a*cos(t);
Y = O_start.a*sin(t);
w = atan2(y2-y1,x2-x1);
x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
y = (y1+y2)/2 + X*sin(w) + Y*cos(w);
plot(x,y,'b', LineWidth=1.2)
hold on 
scatter(-1159.6, 12243.8,'filled')
grid on
axis equal
sist_can(1)
xlim([-5*10^4 2*10^4])
ylim([-2*10^4 2*10^4])
xlabel('e [km]')
ylabel('p [km]')
legend("Orbita iniziale")
hold off

%% Altra rappresentazione
figure
plotOrbit_plane(O_start,0,2*pi,dth,'b');
hold on 
axis equal
scatter(rr_start(1),rr_start(2),'filled')
legend("Orbita iniziale")
xlim([-2*10^4 3*10^4])
ylim([-1*10^4 5*10^4])

%% Rappresentazione grafica orbita finale in ECI
Terra_3D(r_terra)
set(gcf, 'Name', 'Orbita finale', 'NumberTitle', 'off');
scatter3(rr_end(1),rr_end(2),rr_end(3))
hold on
plotOrbit(O_end,0,2*pi,dth,'r');

%% Rappresentazione grafica orbita finale nel sistema perifocale
x1 = rp2;
y1 = 0;
x2 = -ra2;
y2 = 0;
t = linspace(0,2*pi);
X = O_end.a*cos(t);
Y = b2*sin(t);
w = atan2(y2-y1,x2-x1);
x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
y = (y1+y2)/2 + X*sin(w) + Y*cos(w);
plot(x,y,'r', LineWidth=1.2)
hold on 
scatter(rp2, 0,'filled')
grid on
axis equal
sist_can(1)
xlim([-7*10^4 2*10^4])
ylim([-4*10^4 4*10^4])
xlabel('e [km]')
ylabel('p [km]')
legend("Orbita finale")
hold off

%% Altra rappresentazione
figure
plotOrbit_plane(O_end,0,2*pi,dth,'r');
hold on 
scatter(rr_end(1),rr_end(2),'filled')
axis equal
legend("Orbita finale")
xlim([-2*10^4 6*10^4])
ylim([-4*10^4 3*10^4])

%% Confronto le due orbite 
Terra_3D(r_terra)
set(gcf, 'Name', 'Orbita di partenza e orbita di arrivo', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,0,2*pi,dth,'b');
plotOrbit(O_end,0,2*pi,dth,'r');
legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita finale')

%% Strategia 1 
% Trasferimento standard: cambio piano, modifico anomalia pericentro, effettuo il trasferimento bitangente 

% Calcolo i costi

% 1) cambio piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end); 
% 2) modifico anomalia pericentro
[dv2, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);
% 3) modifico orbita con bitangente PERICENTRO-APOCENTRO
[dv3, dv4, dt4, O_bt] = bitangentTransfer(O_cper, O_end, 'pa'); % testando tutte e 4 le possibilità pa è la più economica

% Costo Totale
DV_1=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);


% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino al punto di cambio di pi<ano
dt1 = TOF(O_start, th_start, th_cp); 
% tempo dalla posizione di cambio di piano fino alla posizione di cambio anomalia di pericentro
dt2 = TOF(O_cp, th_cp, th_best(1)); 
% tempo dalla posizione dopo il cambio anomalia di pericentro fino al pericentro:
dt3 = TOF(O_cper,th_best(2),2*pi); 
% tempo della manovra bitangente dt4
% tempo di attesa sull'orbita finale fino al punto d'arrivo
dt5 = TOF(O_end,pi,th_end); 

% Tempo totale 
DT_1= dt1 + dt2 + dt3 + dt4 + dt5;

dtime = seconds(DT_1);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_1);
fprintf("Tempo impiegato: %s \n", dtime);

% rappresentazione grafica
Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 1: cambio piano - cambio pericentro - trasferimento bitangente', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3));
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,th_best(1),dth,'m');
plotOrbit(O_cper,th_best(2),0,dth,'g');
plotOrbit(O_bt,0,pi,dth,'c');
plotOrbit(O_end,pi,th_end,dth,'r');

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita modificata di piano','Orbita modificata anomalia pericentro','orbita bitangente','Orbita finale')

%% Strategia 2 
% Effettuo il trasferimento bitangente su un orbita ausiliaria della stessa forma dell'orbita finale, 
% cambio piano e poi anomalia del pericentro.
% Soluzione che impiega meno tempo 

% Calcolo i costi 

% 1) Trasferimento bitangente su orbita ausiliaria 
[dv1, dv2, dt1, O_bt,th0, thf, O_aus] = bitangentTransfer(O_start, O_end, 'pp');
% 2) cambio piano
[dv3, th_cp, O_cp] = changeOrbitalPlane(O_aus, O_end);
% 3) modifico anomalia pericentro
[dv4, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);

% Costo totale
DV_2=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);

% Calcolo il tempo impiegato

% tempo dalla posizione iniziale fino al pericentro
dt2 = TOF(O_start, th_start, 2*pi); 
% tempo della manovra bitangente dt1
% tempo sull'orbita ausiliaria 
dt3 = TOF(O_aus,0,th_cp); 
% tempo dalla posizione di cambio di piano alla posizione di cambio di anomalia di pericentro
dt4 = TOF(O_cp, th_cp, th_best(1));
% tempo di attesa sull'orbita finale fino al punto d'arrivo (O_cper=O_end)
dt5 = TOF(O_cper,th_best(2),th_end);

% Tempo totale 
DT_2= dt1 + dt2 + dt3 + dt4 + dt5;

dtime = seconds(DT_2);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_2);
fprintf("Tempo impiegato: %s \n", dtime);

% rappresentazione grafica
Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 2: trasferimento bitangente - cambio piano - cambio pericentro', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,2*pi,dth,'b');
plotOrbit(O_bt,0,pi,dth,'c');
plotOrbit(O_aus,0,th_cp,dth,'g');
plotOrbit(O_cp,th_cp,th_best(1),dth,'m');
plotOrbit(O_cper,th_best(2),th_end,dth,'r');

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita bitangente','Orbita ausiliaria', 'Orbita cambio piano', 'Orbita cambio pericentro coincidende con finale')
 
%% Strategia 3  
% Cambio il piano, effettuo il trasferimento bitangente e nel mentre cambio l'anomalia
% del pericentro.

% Calcolo i costi

% 1) cambio piano 
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);  
% 2) trasferimento bitangente 
[dv3, dv4, dt_t, O_bt] = bitangentTransfer(O_cp, O_end, 'aa');
% 3) modifico anomalia pericentro
[dv2, thi, thf, th_best, O_cper] = change_pericenter_arg(O_bt, O_end.om, th_cp);

% Costo totale 
DV_4=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);

% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino al punto di cambio di piano
dt1 = TOF(O_start, th_start, th_cp); 
% tempo dalla posizione di cambio di piano all'apocentro
dt2 = TOF(O_cp, th_cp, pi);
% tempo sulla bitangente ausiliaria
dt3 = TOF(O_bt, 0, th_best(1)-pi); 
% tempo dalla posizione dopo il cambio di anomalia di pericentro fino all'apocentro:
dt4 = TOF(O_cper,th_best(2)-pi,pi); 
% tempo di attesa sull'orbita finale fino al punto d'arrivo
dt5 = TOF(O_end,pi,th_end); 


% Tempo totale
DT_4 = dt1 + dt2 + dt3 + dt4 + dt5;

dtime = seconds(DT_4);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_4);
fprintf("Tempo impiegato: %s \n", dtime);

% rappresentazione grafica
Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 3: cambio piano - trasferimento bitangente - cambio pericentro', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,pi,dth,'m');
plotOrbit(O_bt,0,th_best(1)-pi,dth,'c');
plotOrbit(O_cper,th_best(2)-pi,pi,dth,'g');
plotOrbit(O_end,pi,th_end,dth,'r');

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita cambio piano','Orbita bitangente','Orbita cper','Orbita finale')

%% Strategia 4 
% Durante il trasferimento bitangente cambio il piano e l'anomalia del pericentro.
% Compromesso tra costo e tempo

% Calcolo i costi

% 1) trasferimento bitangente
[dv1, dv2, dt_tot, O_bt, th0, thf, orbit_arr] = bitangentTransfer(O_start, O_end, 'ap');
% 2) modifico piano
[dv3, th_cp, O_cp] = changeOrbitalPlane(O_bt, O_end); 
% 3) cambio anomalia pericentro
[dv4, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);


% Costo totale 
DV_3 = abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);

% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino all'apocentro
dt1 = TOF(O_start, th_start, pi); 
% tempo sulla bitangente ausiliaria
dt2 = TOF(O_bt, pi, th_cp);
% tempo dalla posizione di cambio di piano alla posizione di cambio anomalia di pericentro
dt3 = TOF(O_cp, th_cp, th_best(1)); 
% tempo dalla posizione dopo il cambio di anomalia di pericentro fino al pericentro:
dt4 = TOF(O_cper,th_best(2),2*pi); 
% tempo di attesa sull'orbita finale fino al punto d'arrivo


% Tempo totale
DT_3 = dt1 + dt2 + dt3 + dt4;

dtime = seconds(DT_3);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_3);
fprintf("Tempo impiegato: %s \n", dtime);

% rappresentazione grafica
Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 3: cambio piano e anomalia del pericentro durante il trasferimento bitangente', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start, th_start, pi, dth, 'b');      
plotOrbit(O_bt, pi, th_cp, dth, 'c');
plotOrbit(O_cp, th_cp, th_best(1), dth, 'm'); 
plotOrbit(O_cper, th_best(2), 2*pi,dth,'g');
plotOrbit(O_end, 0, th_end,dth,'r--');

legend('Attrattore','Partenza', 'Arrivo','Orbita iniaziale','Bitangente ausiliaria','Orbita cambio piano','Orbita cambio pericentro', 'Orbita finale')

%% Strategia circolare
% strategia più economica 

% definizione di un'orbita circolare
O_circ = O_start;
O_circ.a = ra1;
O_circ.e = 0;

[dv3, dv4, dt4, O_bt] = bitangentTransfer(O_start, O_circ, 'aa');  % dv4 nullo 
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_bt, O_end); 
[dv2, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp); % dv2 nullo 
[dv33, dv44, dt44, O_btt] = bitangentTransfer(O_cper, O_end, 'pa');

% Costo totale 
DV_1=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+abs(dv33)+abs(dv44);

dt1 = TOF(O_start, th_start, pi); 
dt2 = TOF(O_bt, 0, th_cp); 
dt3 = TOF(O_cp,th_cp,th_best(1)); 
dt_4 = TOF(O_cper,th_best(2),0); 
dt5 = TOF(O_btt,0,pi); 
dt6 = TOF(O_end,pi,th_end); 


% Tempo totale
DT_1= dt1 + dt2 + dt3 + dt_4 + dt5 + dt6;

dtime = seconds(DT_1);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_1);
fprintf("Tempo impiegato: %s \n", dtime);

% rappresentazione grafica
Terra_3D(r_terra)
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3));
plotOrbit(O_start,th_start,pi,dth,'b');
plotOrbit(O_bt,0,th_cp,dth,'c');
plotOrbit(O_cp,th_cp,th_best(1),dth,'k');
plotOrbit(O_cper,th_best(2),0,dth,'g');
plotOrbit(O_btt,0,pi,dth,'m');
plotOrbit(O_end,pi,th_end,dth,'r');

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','orbita trasf 1','Orbita modificata di piano','Orbita modificata anomalia pericentro','orbita trasf 2','Orbita finale')

%% ******************* Trasferimenti biellittici *********************

% Strategia standard applicata al trasferimento biellitico: cambio piano, cambio anomalia del pericentro, trasferimento biellittico.
% Cerco raggio dell'orbita di trasferimento migliore per minimizzare il costo 

[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);
[dv2, thi_cper, thf_cper, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);

dVtot_v = [];
%ra_t_v =[62589:1:200000];
ra_t_v =[40000:1:200000];
n = length(ra_t_v);
for i = 1: n
    ra_t = ra_t_v(i);
    [dv3, dv4, dv5, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_cper,O_end, ra_t); 
    DV_5=abs(dv3)+abs(dv4)+abs(dv5);
    dVtot_v =[ dVtot_v, DV_5];
end
% 
dv_STD_bitan = 0.7304; % Costo della bitangente nella STD
figure
plot(ra_t_v, dVtot_v, 'k', LineWidth=1.3)
hold on
yline(dv_STD_bitan, 'r', LineWidth=1.3)
xlabel('r\_at (km)')
ylabel('\DeltaV (km/s)')
grid on
title('Costo in funzione del raggio di apocentro dell''orbita di trasferimento' )
legend ('Cost della biellittica al variare di r_{at}', 'Costo della Bitangente nella STD')
%% Strategia 5.1

% Calcolo i costi 

% 1) piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end); %la funzione calcola il costo solo in funzione del piano di arrivo, l'orbita che restituisce non è O_end
% 2) modifico anomalia pericentro
[dv2, thi_cper, thf_cper, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);
% 3) trasferimento biellittico  
% ra_t = 10* O_cper.a*(1+O_cper.e);
 ra_t = 64000;

[dv3, dv4, dv5, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_cper,O_end, ra_t); 

% Costo totale
DV_5=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+abs(dv5);


% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino al punto di cambio di piano
dt1 = TOF(O_start, th_start, th_cp); 
% tempo dalla posizione di cambio piano al punto di cambio anomalia di pericentro 
dt2 = TOF(O_cp, th_cp, th_best(1));
% tempo dalla posizione dopo il cambio di anomalia di pericentro fino al pericentro 
dt3 = TOF(O_cper, th_best(2), 0); 
% tempo sulla prima ellisse 
dt4 = TOF(O_biel1,0,pi); 
% tempo sulla seconda ellisse
dt5 = TOF(O_biel2,pi,0); 

% Tempo totale
DT_5 = dt1 + dt2 + dt3 + dt4 + dt5; 

dtime = seconds(DT_5);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_5);
fprintf("Tempo impiegato: %s \n", dtime);


Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 6: CP - CPer - TBT', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b'); 
plotOrbit(O_cp,th_cp,th_best(1),dth,'m'); 
plotOrbit(O_cper,th_best(2),0,dth,'g');
plotOrbit(O_biel1,0,pi,dth,'c');
plotOrbit(O_biel2,pi,0,dth,'k');
plotOrbit(O_end,0,th_end,dth,'r--');
legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita cambio piano','Orbita cambio anomalia pericentro','Orbita ellisse 1','Orbita ellisse 2', 'Orbita finale')

%% Strategia 5.2 - raggio
% Effettuto primo impulso e mi porto sulla prima ellisse da dove cambio il piano, mi porto poi sulla seconda 
% ellise e poi modifico l'anomalia del pericentro

% più economica 

dVtot_v = [];
% ra_t_v =[62589:1:180000];
ra_t_v =[30000:1:180000];
n = length(ra_t_v);
for i = 1: n
    ra_t=ra_t_v(i);
[dv1, dv2, dv3, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_start,O_end, ra_t);
[dv4, th_cp, O_cp] = changeOrbitalPlane(O_biel1, O_end);
O_biel2.i=O_cp.i;
O_biel2.OM=O_cp.OM;
O_biel2.om= O_cp.om;

[dv5, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_biel2, O_end.om, th_cp);

% Costo totale
DV_5=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+ abs(dv5);
dVtot_v =[ dVtot_v, DV_5];

end

figure
plot(ra_t_v, dVtot_v, 'k', LineWidth=1.3)
xlabel('r\_at (km)')
ylabel('\DeltaV (km/s)')
title('Costo in funzione del raggio di apocentro dell''orbita di trasferimento' )

%% Strategia 5.2   
% Cambio il piano nella prima ellisse
% più economica 
 
% Calcolo i costi 

% 1) trasferimento biellittico 
ra_t = 43500;
[dv1, dv2, dv3, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_start,O_end, ra_t);
% 2) cambio piano
[dv4, th_cp, O_cp] = changeOrbitalPlane(O_biel1, O_end);
% 3) finisco biellittico
O_biel2.i=O_cp.i;
O_biel2.OM=O_cp.OM;
O_biel2.om= O_cp.om;
% 4) cambio anomalia pericentro
[dv5, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_biel2, O_end.om, th_cp);

% Costo totale
DV_5=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+ abs(dv5);

 
% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino al pericentro 
dt1 = TOF(O_start, th_start, 2*pi); 
% tempo sulla prima ellisse fino al punto di cambio piano 
dt2 = TOF(O_biel1, 0, th_cp);
% tempo dal punto di cambio piano fino all'apocentro 
dt3 = TOF(O_cp, th_cp, pi); 
% tempo sulla seconda ellisse fino al punto di cambio anomalia del pericentro 
dt4 = TOF(O_biel2,pi,th_best(1)); 
% tempo dal punto di cambio anomalia del pericentro fino al punto finale  
dt5 = TOF(O_cper,th_best(2),th_end); 



% Tempo totale
DT_5 = dt1 + dt2 + dt3 + dt4 + dt5; 

dtime = seconds(DT_5);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_5);
fprintf("Tempo impiegato: %s \n", dtime);

Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 5: TBT - CP - CPer - TBT', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,2*pi,dth,'b');
plotOrbit(O_biel1,0,th_cp,dth, 'm');
plotOrbit(O_cp,th_cp,pi,dth, 'c');
plotOrbit(O_biel2,pi,th_best(1),dth, 'k'); 
plotOrbit(O_cper,th_best(2),th_end, dth, 'g');
plotOrbit(O_end,0,th_end,dth, 'r--');
legend ('Attrattore','Partenza','Arrivo','Orbita iniziale','Prima ellisse','Orbita cambio piano', 'Orbita cambio pericentro', 'Seconda ellisse', 'Orbita finale')

%% Valutazione opzione migliore

% caso 1: cambio piano nella prima ellisse

dVtot_v1 = [];
% ra_t_v =[62589:1:180000];
ra_t_v1 =[300000:1:450000];
n1 = length(ra_t_v1);
for i = 1: n1
    ra_t1=ra_t_v1(i);
[dv11, dv22, dv33, d_t11, d_t22, O_biel11, O_biel22] = biellipticTransfer(O_start,O_end, ra_t1);
[dv44, th_cpp, O_cpp] = changeOrbitalPlane(O_biel11, O_end);
O_biel22.i=O_cpp.i;
O_biel22.OM=O_cpp.OM;
O_biel22.om= O_cpp.om;

[dv55, th_ii, th_ff, th_bestt, O_cperr] = change_pericenter_arg(O_biel22, O_end.om, th_cpp);
DV_66=abs(dv11)+abs(dv22)+abs(dv33)+abs(dv44)+ abs(dv55);

dVtot_v1 =[ dVtot_v1, dv44];

end

% caso 2: cambio piano nella seconda ellisse

dVtot_v = [];
% ra_t_v =[62589:1:180000];
ra_t_v =[300000:1:450000];
n = length(ra_t_v);
for i = 1: n
    ra_t=ra_t_v(i);
[dv1, dv2, dv3, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_start,O_end, ra_t);

[dv4, th_cp, O_cp] = changeOrbitalPlane(O_biel2, O_end);

[dv5, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);

DV_6=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+ abs(dv5);

dVtot_v =[ dVtot_v, dv4];

end



figure
plot(ra_t_v1, dVtot_v1, 'r', LineWidth=1.3)
hold on 
plot(ra_t_v, dVtot_v, 'k', LineWidth=1.3)
xlabel('r\_at (km)')
ylabel('\DeltaV (km/s)')
title('Costo manovara cambio piano in funzione di r\_at' )
legend('Caso 1', 'Caso 2')

%% Strategia 6
% Cambio il piano nella seconda ellisse
% Calcolo i costi 

% 1) trasferimento biellittico 
% ra_t = 2 * O_start.a*(1+O_start.e); 
ra_t = 43500;
[dv1, dv2, dv3, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_start,O_end, ra_t);
% 2 ) cambio piano
[dv4, th_cp, O_cp] = changeOrbitalPlane(O_biel2, O_end);
% 3) cambio anomalia pericentro
[dv5, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);

% Costo totale
DV_6=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+ abs(dv5);

% Calcolo il tempo impiegato 

% ttempo dalla posizione iniziale fino al pericentro 
dt1 = TOF(O_start, th_start, 2*pi); 
% tempo sulla prima ellisse fino all'apocentro 
dt2 = TOF(O_biel1, 0, pi);
% tempo sulla seconda ellisse fino al punto di cambio piano 
dt3 = TOF(O_biel2, pi, th_cp); 
% tempo dal punto di cambio piano fino al punto di cambio di anomalia di pericentro
dt4 = TOF(O_cp,th_cp,th_best(1)); 
% tempo di attesa sull'orbita di cambio di anomalia di pericentro
dt5 = TOF(O_cper,th_best(2),th_end); 


% Tempo totale
DT_6 = dt1 + dt2 + dt3 + dt4 + dt5; 

dtime = seconds(DT_6);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_6);
fprintf("Tempo impiegato: %s \n", dtime);

 

Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 5: TBT - CP - CPer - TBT', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,2*pi,dth,'b');
plotOrbit(O_biel1,2*pi,pi,dth, 'm');
plotOrbit(O_biel2,pi,th_cp,dth, 'k');
plotOrbit(O_cp,th_cp,th_best(1),dth, 'c');
plotOrbit(O_cper,th_best(2),th_end, dth, 'g');
plotOrbit(O_end,0,th_end,dth, 'r--');
legend ('Attrattore','Partenza','Arrivo','Orbita iniziale','Prima ellisse','Seconda ellisse','Orbita cambio piano', 'Orbita cambio pericentro', 'Orbita finale')


