%% Scenario 1
% Trasferimento geocentrico - da GTO ad orbita di parcheggio assegnata
clear all
close all
clc

% Parametro gravitazionale della Terra
mu=398600;  % [km^3/s^2]
% Raggio Terra 
r_terra = 6378.388; % [km]

dth=0.001;    % passo [s]

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

%% Rappresentazione grafica 
Terra_3D(r_terra)
set(gcf, 'Name', 'Orbita di partenza e orbita di arrivo', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,0,2*pi,dth,'b');
plotOrbit(O_end,0,2*pi,dth,'r');
legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita finale')

%% Strategia 1   CP - CPer - TEB   VERSIONE 1 aa  
% QUESTA 
% Strategia standard: cambio piano, modifico anomalia pericentro, effettuo il trasferimento bitangente 

% Calcolo i costi

% 1) cambio piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end); 
% 2) modifico anomalia pericentro
[dv2, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om-pi, th_cp);
% 3) modifico orbita con bitangente 
[dv3, dv4, dt4, O_bt] = bitangentTransfer(O_cper, O_end, 'aa'); % testando tutte e 4 le possibilità conviene aa 
 
% Costo Totale
DV_1=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);



% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino al punto di cambio di piano
dt1 = TOF(O_start, th_start, th_cp); 
% tempo dalla posizione di cambio di piano alla posizione di cambio di anomalia di pericentro
dt2 = TOF(O_cp, th_cp, th_best(1)); 
% tempo dalla posizione dopo il cambio di anomalia di pericentro fino al pericentro:
dt3 = TOF(O_cper,th_best(2),pi); 
% tempo della manovra bitangente dt4
% tempo di attesa sull'orbita finale fino al punto d'arrivo
dt5 = TOF(O_end,pi,th_end); 

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
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,th_best(1),dth,'m'); % Orbita post cambio piano
plotOrbit(O_cper,th_best(2),pi,dth,'g');
plotOrbit(O_bt,0,pi,dth,'c');
plotOrbit(O_end,pi,0,dth,'r');

% plotOrbit(O_start,0,2*pi,dth, ['--','b']);
% h = findobj(gca, 'Type', 'line', '-not', 'HandleVisibility', 'off');
% h(1).HandleVisibility = 'off';  % Nasconde dalla legenda
% plotOrbit(O_end,0,2*pi,dth,['--','r']);
% h = findobj(gca, 'Type', 'line', '-not', 'HandleVisibility', 'off');
% h(1).HandleVisibility = 'off';  % Nasconde dalla legenda

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita modificata di piano','Orbita modificata anomalia pericentro','orbita bitangente','Orbita finale')

%% Strategia 1   CP - CPer - TEB     VERSIONE 2 ap 
% Strategia standard: cambio piano, modifico anomalia pericentro, effettuo il trasferimento bitangente 

% Calcolo i costi

% 1) cambio piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end); 
% 2) modifico anomalia pericentro
[dv2, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);
% 3) modifico orbita con bitangente 
[dv3, dv4, dt4, O_bt] = bitangentTransfer(O_cper, O_end, 'ap'); % testando tutte e 4 le possibilità conviene aa 
 
% Costo Totale
DV_1=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);



% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino al punto di cambio di piano
dt1 = TOF(O_start, th_start, th_cp); 
% tempo dalla posizione di cambio di piano alla posizione di cambio di anomalia di pericentro
dt2 = TOF(O_cp, th_cp, th_best(1)); 
% tempo dalla posizione dopo il cambio di anomalia di pericentro fino al pericentro:
dt3 = TOF(O_cper,th_best(2),pi); 
% tempo della manovra bitangente dt4
% tempo di attesa sull'orbita finale fino al punto d'arrivo
%dt5 = TOF(O_end,pi,th_end); non lo metto perché con l'ultimo impulso sono già lì
DT_1= dt1 + dt2 + dt3 + dt4; %+ dt5;

dtime = seconds(DT_1);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_1);
fprintf("Tempo impiegato: %s \n", dtime);



Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 1: cambio piano - cambio pericentro - trasferimento bitangente', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,th_best(1),dth,'m'); % Orbita post cambio piano
plotOrbit(O_cper,th_best(2),pi,dth,'g');
plotOrbit(O_bt,pi,0,dth,'c');
plotOrbit(O_end,0,2*pi,dth,'r--');


legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita modificata di piano','Orbita modificata anomalia pericentro','orbita bitangente','Orbita finale')


%% Strategia 1   CP - CPer - TEB    VERSIONE 3 pa (più economica) 
% QUESTA
% Strategia standard: cambio piano, modifico anomalia pericentro, effettuo il trasferimento bitangente 

% Calcolo i costi

% 1) cambio piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end); 
% 2) modifico anomalia pericentro
[dv2, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);
% 3) modifico orbita con bitangente 
[dv3, dv4, dt4, O_bt] = bitangentTransfer(O_cper, O_end, 'pa'); % testando tutte e 4 le possibilità conviene aa 
 
% Costo Totale
DV_1=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);



% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino al punto di cambio di piano
dt1 = TOF(O_start, th_start, th_cp); 
% tempo dalla posizione di cambio di piano alla posizione di cambio di anomalia di pericentro
dt2 = TOF(O_cp, th_cp, th_best(1)); 
% tempo dalla posizione dopo il cambio di anomalia di pericentro fino al pericentro:
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
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,th_best(1),dth,'m'); % Orbita post cambio piano
plotOrbit(O_cper,th_best(2),0,dth,'g');
plotOrbit(O_bt,0,pi,dth,'c');
plotOrbit(O_end,pi,0,dth,'r');

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita modificata di piano','Orbita modificata anomalia pericentro','orbita bitangente','Orbita finale')


%% Strategia 1   CP - CPer - TEB    VERSIONE 4 pp buona a livello di tempo,  pessima per costo
% Strategia standard: cambio piano, modifico anomalia pericentro, effettuo il trasferimento bitangente 

% Calcolo i costi

% 1) cambio piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end); 
% 2) modifico anomalia pericentro
[dv2, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om-pi, th_cp);
% 3) modifico orbita con bitangente 
[dv3, dv4, dt4, O_bt] = bitangentTransfer(O_cper, O_end, 'pp'); % testando tutte e 4 le possibilità conviene aa 
 
% Costo Totale
DV_1=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);



% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino al punto di cambio di piano
dt1 = TOF(O_start, th_start, th_cp); 
% tempo dalla posizione di cambio di piano alla posizione di cambio di anomalia di pericentro
dt2 = TOF(O_cp, th_cp, th_best(1)); 
% tempo dalla posizione dopo il cambio di anomalia di pericentro fino al pericentro:
dt3 = TOF(O_cper,th_best(2),2*pi); 
% tempo della manovra bitangente dt4
% tempo di attesa sull'orbita finale fino al punto d'arrivo
% dt5 = TOF(O_end,pi,0); 

% Tempo totale 
DT_1= dt1 + dt2 + dt3 + dt4; %+ dt5;

dtime = seconds(DT_1);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_1);
fprintf("Tempo impiegato: %s \n", dtime);



Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 1: cambio piano - cambio pericentro - trasferimento bitangente', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,th_best(1),dth,'m'); % Orbita post cambio piano
plotOrbit(O_cper,th_best(2),2*pi,dth,'g');
plotOrbit(O_bt,0,pi,dth,'c');
plotOrbit(O_end,0,2*pi,dth,'r--');

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita modificata di piano','Orbita modificata anomalia pericentro','orbita bitangente','Orbita finale')

%% Strategia 2  TEB - CP - CPer VERSIONE 1 aa
% Effettuo Trasferimento bitangente su un orbita ausiliaria, cambio piano e
% poi anomalia del pericentro 

% Calcolo i costi 

% 1) Trasferimento bitangente su orbita ausiliaria 
[dv1, dv2, dt1, O_bt,th0, thf, O_aus] = bitangentTransfer(O_start, O_end, 'aa');
% 2) cambio piano
[dv3, th_cp, O_cp] = changeOrbitalPlane(O_aus, O_end);
% 3) modifico anomalia pericentro
[dv4, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);

% Costo totale
DV_2=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);

% Calcolo il tempo impiegato

% tempo dalla posizione iniziale fino all'apocentro
dt2 = TOF(O_start, th_start, pi); 
% tempo della manovra bitangente dt1
% tempo sull'orbita ausiliaria 
dt3 = TOF(O_aus,pi,th_cp+2*pi); 
% tempo dalla posizione di cambio di piano alla posizione di cambio di anomalia di pericentro
dt4 = TOF(O_cp, th_cp+2*pi, th_best(1));
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
plotOrbit(O_start,th_start,pi,dth,'b');
plotOrbit(O_bt,0,pi,dth,'c');
plotOrbit(O_aus,pi,th_cp+2*pi,dth,'g');
plotOrbit(O_cp,th_cp+2*pi,th_best(1),dth,'m');
plotOrbit(O_cper,th_best(2),2*pi,dth,'r');

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita bitangente','Orbita ausiliaria', 'Orbita cambio piano', 'Orbita cambio pericentro coincidende con finale')
 
%% Strategia 2  TEB - CP - CPer VERSIONE 2 pp
% Effettuo Trasferimento bitangente su un orbita ausiliaria, cambio piano e
% poi anomalia del pericentro 

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

% tempo dalla posizione iniziale fino all'apocentro
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
 
%% Strategia 2  TEB - CP - CPer VERSIONE 3 ap
% Effettuo Trasferimento bitangente su un orbita ausiliaria, cambio piano e
% poi anomalia del pericentro 

% Calcolo i costi 

% 1) Trasferimento bitangente su orbita ausiliaria 
[dv1, dv2, dt1, O_bt,th0, thf, O_aus] = bitangentTransfer(O_start, O_end, 'ap');
% 2) cambio piano
[dv3, th_cp, O_cp] = changeOrbitalPlane(O_aus, O_end);
% 3) modifico anomalia pericentro
[dv4, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);

% Costo totale
DV_2=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);

% Calcolo il tempo impiegato

% tempo dalla posizione iniziale fino all'apocentro
dt2 = TOF(O_start, th_start, pi); 
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

Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 2: trasferimento bitangente - cambio piano - cambio pericentro', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,pi,dth,'b');
plotOrbit(O_bt,pi,0,dth,'c');
plotOrbit(O_aus,0,th_cp,dth,'g');
plotOrbit(O_cp,th_cp,th_best(1),dth,'m');
plotOrbit(O_cper,th_best(2),th_end,dth,'r');

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita bitangente','Orbita ausiliaria', 'Orbita cambio piano', 'Orbita cambio pericentro coincidende con finale')
 
%% Strategia 2  TEB - CP - CPer VERSIONE 4 pa (più economica)
% QUESTA
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
plotOrbit(O_start,th_start,2*pi,dth,'b');
plotOrbit(O_bt,0,pi,dth,'c');
plotOrbit(O_aus,pi,th_cp,dth,'g');
plotOrbit(O_cp,th_cp,th_best(1),dth,'m');
plotOrbit(O_cper,th_best(2),th_end,dth,'r');


legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita bitangente','Orbita ausiliaria', 'Orbita cambio piano', 'Orbita cambio pericentro coincidende con finale')
 
%% Strategia 3 TBE - CP - CPer VERSIONE 1 aa  (più economica)
% QUESTA
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
plotOrbit(O_start, th_start, pi, dth, 'b');          
plotOrbit(O_bt, 0, th_cp, dth, 'c');
plotOrbit(O_cp, th_cp, th_best(1), dth, 'm'); 
plotOrbit(O_cper, th_best(2), pi,dth,'g');
plotOrbit(O_end, pi, th_end,dth,'r');

legend('Attrattore','Partenza', 'Arrivo','Orbita iniaziale','Bitangente ausiliaria','Orbita cambio piano','Orbita cambio pericentro', 'Orbita finale')

%% Strategia 3 TBE - CP - CPer VERSIONE 2 pp
% Durante il trasferimento bitangente cambio il piano e l'anomalia del pericentro 

% Calcolo i costi

% 1) trasferimento bitangente
[dv1, dv2, dt_tot, O_bt, th0, thf, orbit_arr] = bitangentTransfer(O_start, O_end, 'pp');
% 2) modifico piano
[dv3, th_cp, O_cp] = changeOrbitalPlane(O_bt, O_end); 
% 3) cambio anomalia pericentro
[dv4, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om-pi, th_cp);


% Costo totale 
DV_3 = abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);

% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino all'apocentro
dt1 = TOF(O_start, th_start, 2*pi); 
% tempo sulla bitangente ausiliaria
dt2 = TOF(O_bt, 0, th_cp);
% tempo dalla posizione di cambio di piano alla posizione di cambio anomalia di pericentro
dt3 = TOF(O_cp, th_cp, th_best(1)); 
% tempo dalla posizione dopo il cambio di anomalia di pericentro fino al pericentro:
dt4 = TOF(O_cper,th_best(2),pi); 
% tempo di attesa sull'orbita finale fino al punto d'arrivo
%dt5 = TOF(O_end,pi,th_end); 

% Tempo totale
DT_3 = dt1 + dt2 + dt3 + dt4; %+ dt5;

dtime = seconds(DT_3);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_3);
fprintf("Tempo impiegato: %s \n", dtime);

Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 3: cambio piano e anomalia del pericentro durante il trasferimento bitangente', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start, th_start, 2*pi, dth, 'b');    
plotOrbit(O_bt, 0, th_cp, dth, 'c');
plotOrbit(O_cp, th_cp, th_best(1), dth, 'm'); 
plotOrbit(O_cper, th_best(2), pi,dth,'g');
plotOrbit(O_end, 0, th_end,dth,'r--');

legend('Attrattore','Partenza', 'Arrivo','Orbita iniaziale','Bitangente ausiliaria','Orbita cambio piano','Orbita cambio pericentro', 'Orbita finale')

%% Strategia 3 TBE - CP - CPer VERSIONE 3 pa
% Durante il trasferimento bitangente cambio il piano e l'anomalia del pericentro 

% Calcolo i costi

% 1) trasferimento bitangente
[dv1, dv2, dt_tot, O_bt, th0, thf, orbit_arr] = bitangentTransfer(O_start, O_end, 'pa');
% 2) modifico piano
[dv3, th_cp, O_cp] = changeOrbitalPlane(O_bt, O_end); 
% 3) cambio anomalia pericentro
[dv4, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);


% Costo totale 
DV_3 = abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);

% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino all'apocentro
dt1 = TOF(O_start, th_start, 2*pi); 
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
plotOrbit(O_start, th_start, 2*pi, dth, 'b');   
plotOrbit(O_bt, 0, th_cp, dth, 'c');
plotOrbit(O_cp, th_cp, th_best(1), dth, 'm'); 
plotOrbit(O_cper, th_best(2), pi,dth,'g');
plotOrbit(O_end, pi, th_end,dth,'r');

legend('Attrattore','Partenza', 'Arrivo','Orbita iniaziale','Bitangente ausiliaria','Orbita cambio piano','Orbita cambio pericentro', 'Orbita finale')

%% Strategia 3 TBE - CP - CPer VERSIONE 4 ap
% Durante il trasferimento bitangente cambio il piano e l'anomalia del pericentro 

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
%dt5 = TOF(O_end,pi,th_end); 

% Tempo totale
DT_3 = dt1 + dt2 + dt3 + dt4; %+ dt5;

dtime = seconds(DT_3);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_3);
fprintf("Tempo impiegato: %s \n", dtime);

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



%% Strategia 4   CP - TEB - CPer  VERSIONE 1 aa  (più economica)
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
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,pi,dth,'m');
plotOrbit(O_bt,0,th_best(1)-pi,dth,'c');
plotOrbit(O_cper,th_best(2)-pi,pi,dth,'g');
plotOrbit(O_end,pi,th_end,dth,'r');

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita cambio piano','Orbita bitangente','Orbita cper','Orbita finale')

%% Strategia 4   CP - TEB - CPer  VERSIONE 2 pp
% Cambio il piano, effettuo il trasferimento bitangente, cambio l'anomalia
% del pericentro 

% Calcolo i costi

% 1) cambio piano 
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);  
% 2) trasferimento bitangente 
[dv3, dv4, dt_t, O_bt] = bitangentTransfer(O_cp, O_end, 'pp');
% 3) modifico anomalia pericentro
[dv2, thi, thf, th_best, O_cper] = change_pericenter_arg(O_bt, O_end.om-pi, th_cp);

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
set(gcf, 'Name', 'Strategia 3: cambio piano - trasferimento bitangente - cambio pericentro', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,0,dth,'m');
plotOrbit(O_bt,0,th_best(1)-pi,dth,'c');
plotOrbit(O_cper,th_best(2)-pi,pi,dth,'g');
%plotOrbit(O_cper,0,2*pi,dth,'g');
plotOrbit(O_end,0,th_end,dth,'r--');

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita cambio piano','Orbita bitangente','Orbita cper','Orbita finale')

%% Strategia 4   CP - TEB - CPer  VERSIONE 3 pa
% Cambio il piano, effettuo il trasferimento bitangente, cambio l'anomalia
% del pericentro 

% Calcolo i costi

% 1) cambio piano 
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);  
% 2) trasferimento bitangente 
[dv3, dv4, dt_t, O_bt] = bitangentTransfer(O_cp, O_end, 'pa');
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
set(gcf, 'Name', 'Strategia 3: cambio piano - trasferimento bitangente - cambio pericentro', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,0,dth,'m');
plotOrbit(O_bt,0,th_best(1)-pi,dth,'c');
plotOrbit(O_cper,th_best(2)-pi,pi,dth,'g');
plotOrbit(O_end,pi,th_end,dth,'r');

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita cambio piano','Orbita bitangente','Orbita cper','Orbita finale')

%% Strategia 4   CP - TEB - CPer  VERSIONE 4 ap
% Cambio il piano, effettuo il trasferimento bitangente, cambio l'anomalia
% del pericentro 

% Calcolo i costi

% 1) cambio piano 
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);  
% 2) trasferimento bitangente 
[dv3, dv4, dt_t, O_bt] = bitangentTransfer(O_cp, O_end, 'ap');
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
set(gcf, 'Name', 'Strategia 3: cambio piano - trasferimento bitangente - cambio pericentro', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,pi,dth,'m');
plotOrbit(O_bt,pi,th_best(1)-pi,dth,'c');
plotOrbit(O_cper,th_best(2)-pi,0,dth,'g');
plotOrbit(O_end,0,th_end,dth,'r--');

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita cambio piano','Orbita bitangente','Orbita cper','Orbita finale')

%%  Strategia 5  TBT - CP - CPer - TBT VERSIONE 1

% Calcolo i costi 

% 1) trasferimento biellittico 
ra_t = 50 * O_start.a*(1+O_start.e); % n volta l'apocentro di O_start
[dv1, dv2, dv3, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_start,O_end, ra_t);%modifica a,e
% 2) cambio piano
[dv4, th_cp, O_cp] = changeOrbitalPlane(O_biel1, O_end);
% 3) cambio anomalia pericentro
[dv5, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);
% 4) finisco biellittico
O_biel2.i=O_cper.i;
O_biel2.OM=O_cper.OM;
O_biel2.om= O_cper.om;

% Costo totale
DV_5=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+ abs(dv5);

% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino al punto di cambio di piano
dt1 = TOF(O_start, th_start, 2*pi); 
% tempo dalla posizione di cambio di piano all'apocentro
dt2 = TOF(O_biel1, 0, th_cp);
% tempo sulla bitangente ausiliaria
dt3 = TOF(O_cp, th_cp, th_best(1)); 
% tempo dalla posizione dopo il cambio di anomalia di pericentro fino all'apocentro:
dt4 = TOF(O_cper,th_best(2),pi); 
% tempo sulla seconda ellisse
dt5 = TOF(O_biel2,pi,th_end); 
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
plotOrbit(O_start,th_start,2*pi,dth,'b');
plotOrbit(O_biel1,0,th_cp,dth, 'm');
plotOrbit(O_cp,th_cp,th_best(1),dth, 'c');
plotOrbit(O_cper,th_best(2),pi,dth, 'k');
plotOrbit(O_biel2,pi,th_end, dth, 'g');
plotOrbit(O_end,0,th_end,dth, 'r--');
legend ('Attrattore','Partenza','Arrivo','Orbita iniziale','Prima ellisse','Orbita cambio piano', 'Orbita cambio pericentro', 'Seconda ellisse', 'Orbita finale')

%% Strategia 5   TBT - CP - TBT - CPer  VERSIONE 2  (più economica)
 
% Calcolo i costi 

% 1) trasferimento biellittico 
ra_t = 50 * O_start.a*(1+O_start.e); % n volte l'apocentro di O_start
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
plotOrbit(O_start,th_start,2*pi,dth,'b');
plotOrbit(O_biel1,0,th_cp,dth, 'm');
plotOrbit(O_cp,th_cp,pi,dth, 'c');
plotOrbit(O_biel2,pi,th_best(1),dth, 'k'); % qui boh perché la faccio concludere in th_best(1) che non è l'apocentro poi la ruoto e il punto di intersazione con end è proprio il punto di arrivo che coincide con il pericentro e quindi faccio qui ultimo impulso
plotOrbit(O_cper,th_best(2),th_end, dth, 'g');
plotOrbit(O_end,0,th_end,dth, 'r--');
legend ('Attrattore','Partenza','Arrivo','Orbita iniziale','Prima ellisse','Orbita cambio piano', 'Orbita cambio pericentro', 'Seconda ellisse', 'Orbita finale')

%% Strategia 5   TBT - CP - TBT - CPer  VERSIONE 3
 
% Calcolo i costi 

% 1) trasferimento biellittico 
ra_t = 50 * O_start.a*(1+O_start.e); % n volte l'apocentro di O_start
[dv1, dv2, dv3, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_start,O_end, ra_t);%modifica a,e
%2 ) cambio piano
[dv4, th_cp, O_cp] = changeOrbitalPlane(O_biel2, O_end);
% 3) cambio anomalia pericentro
[dv5, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);

% Costo totale
DV_5=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+ abs(dv5);

% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino al punto di cambio di piano
dt1 = TOF(O_start, th_start, 2*pi); 
% tempo dalla posizione di cambio di piano all'apocentro
dt2 = TOF(O_biel1, 0, pi);
% tempo sulla bitangente ausiliaria
dt3 = TOF(O_biel2, pi, th_cp); 
% tempo dalla posizione dopo il cambio di anomalia di pericentro fino all'apocentro:
dt4 = TOF(O_cp,th_cp,th_best(1)); 
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
plotOrbit(O_start,th_start,2*pi,dth,'b');
plotOrbit(O_biel1,2*pi,pi,dth, 'm');
plotOrbit(O_biel2,pi,th_cp,dth, 'k');
plotOrbit(O_cp,th_cp,th_best(1),dth, 'c');
plotOrbit(O_cper,th_best(2),th_end, dth, 'g');
plotOrbit(O_end,0,th_end,dth, 'r--');
legend ('Attrattore','Partenza','Arrivo','Orbita iniziale','Prima ellisse','Seconda ellisse','Orbita cambio piano', 'Orbita cambio pericentro', 'Orbita finale')

%% Strategia    CP - TEB - CPer  
%% scritte sopra a strateria 4
% 1) piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);  
dt1 = TOF(O_start, th_start, th_cp);

[dv3, dv4, dt3, O_bt] = bitangentTransfer(O_cp, O_end, 'pa');

%2)modifico anomalia pericentro
[dv2, thi_cper, thf_cper, th_best, O_cper] = change_pericenter_arg(O_bt, O_end.om, th_cp);
dt2 = TOF(O_cp, th_cp, th_best);



DV_4=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)

Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 5: TBT - CP - CPer - TBT', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,0,dth,'r');
plotOrbit(O_bt,0,th_best(1),dth,'g');
plotOrbit(O_cper,th_best(2),pi,dth,'c');
plotOrbit(O_end,0,th_end,dth,'m');
xlim([-1e5,1e5])
ylim([-1e5,1e5])
zlim([-1e5,1e5])
legend('attrattore','partenza','arrivo','sitema riferimento','orbita di inizio','orbita modificata di piano','orbita bitangente','orbita cp','orbita finale')

%% soluzione   CP - TEB - CPer  
% 1) piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);  
dt1 = TOF(O_start, th_start, th_cp);

[dv3, dv4, dt3, O_bt] = bitangentTransfer(O_cp, O_end, 'ap');

%2)modifico anomalia pericentro
[dv2, thi_cper, thf_cper, th_best, O_cper] = change_pericentre_arg(O_bt, O_end.om, th_cp);
dt2 = TOF(O_cp, th_cp, th_best);



DV_4=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)

sol_1=figure;
sol_1.Name="soluzione CP - TEB - CPer ";
scatter3(0,0,0)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,pi,dth,'r');
plotOrbit(O_bt,pi,th_best(1),dth,'g');
plotOrbit(O_cper,th_best(2),pi+pi,dth,'c');
plotOrbit(O_end,0,th_end,dth,'m');
xlim([-1e5,1e5])
ylim([-1e5,1e5])
zlim([-1e5,1e5])
legend('attrattore','partenza','arrivo','sitema riferimento','orbita di inizio','orbita modificata di piano','orbita bitangente','orbita cp','orbita finale')

%% soluzione   CP - TEB - CPer  
% disegno brutto
% 1) piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);  
dt1 = TOF(O_start, th_start, th_cp);

[dv3, dv4, dt3, O_bt] = bitangentTransfer(O_cp, O_end, 'pp');

%2)modifico anomalia pericentro
[dv2, thi_cper, thf_cper, th_best, O_cper] = change_pericentre_arg(O_bt, O_end.om, th_cp);
dt2 = TOF(O_cp, th_cp, th_best);



DV_4=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)

sol_1=figure;
sol_1.Name="soluzione CP - TEB - CPer ";
scatter3(0,0,0)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,0,dth,'r');
plotOrbit(O_bt,pi + pi,th_best(1),dth,'g');
plotOrbit(O_cper,th_best(2),pi,dth,'c');
plotOrbit(O_end,0,th_end,dth,'m');
xlim([-1e5,1e5])
ylim([-1e5,1e5])
zlim([-1e5,1e5])
legend('attrattore','partenza','arrivo','sitema riferimento','orbita di inizio','orbita modificata di piano','orbita bitangente','orbita cp','orbita finale')

%% Strategia 6   CP - CPer - TBT

% Calcolo i costi 

% 1) piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end); %la funzione calcola il costo solo in funzione del piano di arrivo, l'orbita che restituisce non è O_end
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
plotOrbit(O_start,th_start,th_cp,dth,'b'); % Orbita di partenza
plotOrbit(O_cp,th_cp,th_best(1),dth,'m'); % Orbita post cambio piano
plotOrbit(O_cper,th_best(2),0,dth,'g');
plotOrbit(O_biel1,0,pi,dth,'c');
plotOrbit(O_biel2,pi,0,dth,'k');
plotOrbit(O_end,0,th_end,dth,'r--');
legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita cambio piano','Orbita cambio anomalia pericentro','Orbita ellisse 1','Orbita ellisse 2')

%% Strategia 7    CP - TBT - CPer  VERSIONE 1 (più economica)

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
set(gcf, 'Name', 'Strategia 6: CP - TBT - CPer', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,0,dth, 'm');
plotOrbit(O_biel1,0,pi,dth, 'c');
plotOrbit(O_biel2,pi,th_best(1),dth, 'k');
plotOrbit(O_cper,th_best(2),th_end, dth, 'g');
plotOrbit(O_end,0,th_end,dth, 'r--');
legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita cambio piano','Orbita ellisse 1','Orbita ellisse 2','Orbita cambio anomalia pericentro')

%% Strategia 7    CP - TBT - CPer - TBT VERSIONE 2    

% Calcolo i costi 

% 1) piano orbita
[dv4, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);
% 2)Trasferimento biellittico
ra_t = 2 * O_start.a*(1+O_start.e); % n volta l'apocentro di O_start
[dv1, dv2, dv3, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_cp,O_end, ra_t);
% 3) cambio anomalia del pericentro 
[dv5, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_biel1, O_end.om, th_cp);
% 4) finisco biellittica 
O_biel2.i=O_cper.i;
O_biel2.OM=O_cper.OM;
O_biel2.om= O_cper.om;

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
set(gcf, 'Name', 'Strategia 6: CP - TBT - CPer', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,0,dth, 'm');
plotOrbit(O_biel1,0,th_best(1),dth, 'c'); 
plotOrbit(O_cper,th_best(2),pi, dth, 'g');
plotOrbit(O_biel2,pi,th_end,dth, 'k');
plotOrbit(O_end,0,th_end,dth, 'r--');
legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita cambio piano','Orbita ellisse 1','Orbita cambio anomalia pericentro','Orbita ellisse 2')

%% Strategia 7    CP - TBT - CPer - TBT VERSIONE 3    

% Calcolo i costi 

% 1) piano orbita
[dv4, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);
% 2)Trasferimento biellittico
ra_t = 2 * O_start.a*(1+O_start.e); % n volta l'apocentro di O_start
[dv1, dv2, dv3, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_cp,O_end, ra_t);
% 3) cambio anomalia del pericentro 
[dv5, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_biel1, O_end.om, th_cp);
% 4) finisco biellittica 
O_biel2.i=O_cper.i;
O_biel2.OM=O_cper.OM;
O_biel2.om= O_cper.om;

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
dt5 = TOF(O_cper,th_best(2),0); 
% tempo di attesa sull'orbita finale fino al punto d'arrivo
%dt6 = TOF(O_end,pi,th_end); 


% Tempo totale
DT_7 = dt1 + dt2 + dt3 + dt4 + dt5; %+ dt6;

dtime = seconds(DT_7);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra: %d \n", DV_7);
fprintf("Tempo impiegato: %s \n", dtime);


Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 6: CP - TBT - CPer', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b');
plotOrbit(O_cp,th_cp,0,dth, 'm');
plotOrbit(O_biel1,0,th_i(1),dth, 'c'); 
plotOrbit(O_cper,th_f(1),pi, dth, 'g');
plotOrbit(O_biel2,pi,th_end,dth, 'k');
plotOrbit(O_end,0,th_end,dth, 'r--');
legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita cambio piano','Orbita ellisse 1','Orbita cambio anomalia pericentro','Orbita ellisse 2')
%% 
% individuiamo l'orbita di trasferimento

% e_t = (p2-(1+e2*cos(th_PC))*rp1)/(rp1*(1+e2*cos(th_PC))-p2*cos(th_PC));
% a_t = rp1/(1-e_t);
% p_t = a_t * (1 - e_t^2);
% kepEl_t = [a_t e_t i1 OM1 om1 2*pi]';

th_cp = 3.4563;
p1 = O_start.a*(1-O_start.e^2);
p2 = O_end.a*(1-O_end.e^2);
rp1 = O_start.a*(1-O_start.e);
e_t = (p2-(1+O_end.e*cos(th_cp))*rp1)/(rp1*(1+O_end.e*cos(th_cp))-p2*cos(th_cp));
a_t = rp1/(1-e_t);
p_t = a_t * (1 - e_t^2);
O_tr.a = a_t;
O_tr.e = e_t;
O_tr.OM = 1.514000;
O_tr.om = 3.107000;
O_tr.i = 0.1047;
O_tr.mu = mu; 
 
[dv3, th_cpp, O_cp] = changeOrbitalPlane(O_tr, O_end);
[dv4, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);
[deltav,deltav_r,alpha]=calcolo_velocita_man_combo(O_tr,O_end,th_cp,th_cp);
deltav1_t = abs(sqrt(mu/p_t) * (1 + e_t) - sqrt(mu/p1) * (1 + O_start.e));

% Costo totale
DV_2= abs(dv4) + deltav1_t + abs(deltav);

Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 2: trasferimento bitangente - cambio piano - cambio pericentro', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,0,dth,'b');
plotOrbit(O_tr,0,th_best(1),dth,'c');
plotOrbit(O_cper,th_best(1),0, dth, 'g');
plotOrbit(O_end,0,2*pi,dth,'r');
% plotOrbit(O_cper,th_best(2),th_end,dth,'r');