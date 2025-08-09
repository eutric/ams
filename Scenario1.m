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

%% Rappresentazione grafica orbita iniziale in ECI
Terra_3D(r_terra)
set(gcf, 'Name', 'Orbita di partenza', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
plotOrbit(O_start,0,2*pi,dth,'b')

%% Orbita iniziale nel piano orbitale
% DA AGGIUSTARE
figure
plotOrbit_plane(O_start,0,2*pi,dth,'b')
hold on 
scatter(rr_start(1),rr_start(2),'filled')

%% Rappresentazione grafica orbita finale in ECI
Terra_3D(r_terra)
set(gcf, 'Name', 'Orbita finale', 'NumberTitle', 'off');
scatter3(rr_end(1),rr_end(2),rr_end(3))
hold on
plotOrbit(O_end,0,2*pi,dth,'r')

%% Orbita finale nel piano orbitale
% DA AGGIUSTARE
figure
plotOrbit_plane(O_end,0,2*pi,dth,'r')
hold on 
scatter(rr_end(1),rr_end(2),'filled')

%% Rappresentazione grafica delle due orbite
Terra_3D(r_terra)
set(gcf, 'Name', 'Orbita di partenza e orbita di arrivo', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,0,2*pi,dth,'b')
plotOrbit(O_end,0,2*pi,dth,'r')
legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita finale')

%% Strategia 1 CP - CPer - TEB 
% Trasferimento standard: cambio piano, modifico anomalia pericentro, effettuo il trasferimento bitangente 

% Calcolo i costi

% 1) cambio piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end); 
% 2) modifico anomalia pericentro
[dv2, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);
% 3) modifico orbita con bitangente PERICENTRO-APOCENTRO
[dv3, dv4, dt4, O_bt] = bitangentTransfer(O_cper, O_end, 'pa'); % testando tutte e 4 le possibilità conviene pa 
 
% Costo Totale
DV_1=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);


% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino al punto di cambio di piano
dt1 = TOF(O_start, th_start, th_cp); 
% tempo dalla posizione di cambio di piano fino alla posizione di cambio anomalia di pericentro
dt2 = TOF(O_cp, th_cp, th_best(1)); 
% tempo dalla posizione dopo il cambio anomalia di pericentro fino al pericentro:
dt3 = TOF(O_cper,th_best(2),2*pi); 
% tempo della manovra bitangente dt4
% tempo di attesa sull'orbita finale fino al punto d'arrivo
dt5 = TOF(O_end,pi,0); 

% Tempo totale 
DT_1= dt1 + dt2 + dt3 + dt4 + dt5;

dtime = seconds(DT_1);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_1);
fprintf("Tempo impiegato: %s \n", dtime);



Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 1: cambio piano - cambio pericentro - trasferimento bitangente', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b')
plotOrbit(O_cp,th_cp,th_best(1),dth,'m')
plotOrbit(O_cper,th_best(2),0,dth,'g')
plotOrbit(O_bt,0,pi,dth,'c')
plotOrbit(O_end,pi,0,dth,'r')

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita modificata di piano','Orbita modificata anomalia pericentro','orbita bitangente','Orbita finale')

%% Strategia 2  TEB - CP - CPer 
% Effettuo Trasferimento bitangente su un orbita ausiliaria, cambio piano e
% poi anomalia del pericentro 

% Calcolo i costi 

% 1) Trasferimento bitangente su orbita ausiliaria 
[dv1, dv2, dt1, O_bt,th0, thf, O_aus] = bitangentTransfer(O_start, O_end, 'pa');
% 2) cambio piano
[dv3, th_cp, O_cp] = changeOrbitalPlane(O_aus, O_end);
% 3) modifico anomalia pericentro
[dv4, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);

% Costo totale
DV_2=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);

% Calcolo il tempo impiegato

% tempo dalla posizione iniziale fino all'apocentro
dt2 = TOF(O_start, th_start, 2*pi); 
% tempo della manovra bitangente dt1
% tempo sull'orbita ausiliaria 
dt3 = TOF(O_aus,pi,th_cp); 
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

Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 2: trasferimento bitangente - cambio piano - cambio pericentro', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,2*pi,dth,'b')
plotOrbit(O_bt,0,pi,dth,'c')
plotOrbit(O_aus,pi,th_cp,dth,'g')
plotOrbit(O_cp,th_cp,th_best(1),dth,'m')
plotOrbit(O_cper,th_best(2),th_end,dth,'r')


legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita bitangente','Orbita ausiliaria', 'Orbita cambio piano', 'Orbita cambio pericentro coincidende con finale')

%% Strategia 3 TBE - CP - CPer 
% Durante il trasferimento bitangente cambio il piano e l'anomalia del pericentro 

% Calcolo i costi

% 1) trasferimento bitangente
[dv1, dv2, dt_tot, O_bt, th0, thf, orbit_arr] = bitangentTransfer(O_start, O_end, 'aa');
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
dt2 = TOF(O_bt, 0, th_cp);
% tempo dalla posizione di cambio di piano alla posizione di cambio anomalia di pericentro
dt3 = TOF(O_cp, th_cp, th_best(1)); 
% tempo dalla posizione dopo il cambio di anomalia di pericentro fino al pericentro:
dt4 = TOF(O_cper,th_best(2),pi); 
% tempo di attesa sull'orbita finale fino al punto d'arrivo
dt5 = TOF(O_end,pi,th_end); 

% Tempo totale
DT_3 = dt1 + dt2 + dt3 + dt4 + dt5;

dtime = seconds(DT_3);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_3);
fprintf("Tempo impiegato: %s \n", dtime);

Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 3: cambio piano e anomalia del pericentro durante il trasferimento bitangente', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start, th_start, pi, dth, 'b')           
plotOrbit(O_bt, 0, th_cp, dth, 'c')      
plotOrbit(O_cp, th_cp, th_best(1), dth, 'm') 
plotOrbit(O_cper, th_best(2), pi,dth,'g') 
plotOrbit(O_end, pi, th_end,dth,'r') 

legend('Attrattore','Partenza', 'Arrivo','Orbita iniaziale','Bitangente ausiliaria','Orbita cambio piano','Orbita cambio pericentro', 'Orbita finale')

%% Strategia 4   CP - TEB - CPer  (più economica)
% Cambio il piano, effettuo il trasferimento bitangente, cambio l'anomalia
% del pericentro 

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


Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 4: cambio piano - trasferimento bitangente - cambio pericentro', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b') 
plotOrbit(O_cp,th_cp,pi,dth,'m')
plotOrbit(O_bt,0,th_best(1)-pi,dth,'c')
plotOrbit(O_cper,th_best(2)-pi,pi,dth,'g')
plotOrbit(O_end,pi,th_end,dth,'r')

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita cambio piano','Orbita bitangente','Orbita cper','Orbita finale')


%% Strategia 5   TBT - CP - TBT - CPer 
% trasferimento sulla prima ellisse, cambio piano, cambio anomalia del
% pericentro, trasferimento sulla seconda ellisse e arrivo al punto finale 
 
% Calcolo i costi 

% 1) trasferimento biellittico 
ra_t = 2 * O_start.a*(1+O_start.e); % n volte l'apocentro di O_start
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

% tempo dalla posizione iniziale fino al punto di cambio di piano
dt1 = TOF(O_start, th_start, 2*pi); 
% tempo dalla posizione di cambio di piano all'apocentro
dt2 = TOF(O_biel1, 0, th_cp);
% tempo sulla bitangente ausiliaria
dt3 = TOF(O_cp, th_cp, pi); 
% tempo dalla posizione dopo il cambio di anomalia di pericentro fino all'apocentro:
dt4 = TOF(O_biel2,pi,th_best(1)); 
% tempo sulla seconda ellisse
dt5 = TOF(O_cper,th_best(2),th_end); 
% tempo di attesa sull'orbita finale fino al punto d'arrivo
%dt6 = TOF(O_end,pi,th_end); 


% Tempo totale
DT_5 = dt1 + dt2 + dt3 + dt4 + dt5; %+ dt6;

dtime = seconds(DT_5);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_5);
fprintf("Tempo impiegato: %s \n", dtime);

Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 5: TBT - CP - CPer - TBT', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,2*pi,dth,'b')
plotOrbit(O_biel1,0,th_cp,dth, 'm') 
plotOrbit(O_cp,th_cp,pi,dth, 'c')
plotOrbit(O_biel2,pi,th_best(1),dth, 'k') 
plotOrbit(O_cper,th_best(2),th_end, dth, 'g')
plotOrbit(O_end,0,th_end,dth, 'r--')
legend ('Attrattore','Partenza','Arrivo','Orbita iniziale','Prima ellisse','Orbita cambio piano', 'Orbita cambio pericentro', 'Seconda ellisse', 'Orbita finale')

%% Strategia 6   CP - CPer - TBT
% cambio piano, cambio anomalia del pericentro, trasferimento biellittico 

% Calcolo i costi 

% 1) piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end); 
% 2) modifico anomalia pericentro
[dv2, thi_cper, thf_cper, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);
% 3) trasferimento biellittico  
ra_t = 2 * O_cper.a*(1+O_cper.e);
[dv3, dv4, dv5, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_cper,O_end, ra_t); 

% Costo totale
DV_6=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+abs(dv5);


% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino al punto di cambio di piano
dt1 = TOF(O_start, th_start, th_cp); 
% tempo dalla posizione di cambio di piano all'apocentro
dt2 = TOF(O_cp, th_cp, th_best(1));
% tempo sulla bitangente ausiliaria
dt3 = TOF(O_cper, th_best(2), 0); 
% tempo dalla posizione dopo il cambio di anomalia di pericentro fino all'apocentro:
dt4 = TOF(O_biel1,0,pi); 
% tempo sulla seconda ellisse
dt5 = TOF(O_biel2,pi,0); 
% tempo di attesa sull'orbita finale fino al punto d'arrivo
%dt6 = TOF(O_end,pi,th_end); 


% Tempo totale
DT_6 = dt1 + dt2 + dt3 + dt4 + dt5; %+ dt6;

dtime = seconds(DT_6);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_6);
fprintf("Tempo impiegato: %s \n", dtime);


Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 6: CP - CPer - TBT', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b') 
plotOrbit(O_cp,th_cp,th_best(1),dth,'m') 
plotOrbit(O_cper,th_best(2),0,dth,'g')
plotOrbit(O_biel1,0,pi,dth,'c')
plotOrbit(O_biel2,pi,0,dth,'k')
plotOrbit(O_end,0,th_end,dth,'r--')
legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita cambio piano','Orbita cambio anomalia pericentro','Orbita ellisse 1','Orbita ellisse 2')

%% Strategia 7    CP - TBT - CPer 
% cambio piano, trasferimento biellittico, cambio anomalia del pericentro

% Calcolo i costi 

% 1) piano orbita
[dv4, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);
% 2)Trasferimento biellittico
ra_t = 2 * O_start.a*(1+O_start.e); % n volte l'apocentro di O_start
[dv1, dv2, dv3, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_cp,O_end, ra_t);
% 3) cambio anomalia del pericentro 
[dv5, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_biel2, O_end.om, th_cp);

% Costo totale
DV_7=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+ abs(dv5);
 

% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino al punto di cambio di piano
dt1 = TOF(O_start, th_start, th_cp); 
% tempo dalla posizione di cambio di piano all'apocentro
dt2 = TOF(O_cp, th_cp, 0);
% tempo sulla bitangente ausiliaria
dt3= TOF(O_biel1,0,pi); 
% tempo dalla posizione dopo il cambio di anomalia di pericentro fino all'apocentro:
dt4 = TOF(O_biel2,0,th_best(1)); 
% tempo sulla seconda ellisse
dt5 = TOF(O_cper,th_best(2),th_end); 
% tempo di attesa sull'orbita finale fino al punto d'arrivo
%dt6 = TOF(O_end,pi,th_end); 


% Tempo totale
DT_7 = dt1 + dt2 + dt3 + dt4 + dt5; %+ dt6;

dtime = seconds(DT_7);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_7);
fprintf("Tempo impiegato: %s \n", dtime);


Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 7: CP - TBT - CPer', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b')
plotOrbit(O_cp,th_cp,0,dth, 'm')
plotOrbit(O_biel1,0,pi,dth, 'c') 
plotOrbit(O_biel2,pi,th_best(1),dth, 'k')
plotOrbit(O_cper,th_best(2),th_end, dth, 'g')
plotOrbit(O_end,0,th_end,dth, 'r--')
legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita cambio piano','Orbita ellisse 1','Orbita ellisse 2','Orbita cambio anomalia pericentro')
