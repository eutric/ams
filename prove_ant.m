%% Scenario 1
% Trasferimento geocentrico - da GTO ad orbita di parcheggio assegnata
clear all
close all
clc

% Parametro gravitazionale della Terra
mu=398600;  % [km^3/s^2]
% Raggio Terra 
r_terra = 6378.388; % [km]

dth=0.001;  % passo 

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

%% Strategia 1   CP - CPer - TEB

% Calcolo il costo

% 1) cambio piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end); 


% 2) modifico anomalia pericentro
[dv2, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om-pi, th_cp);

%3)modifico orbita con bitangente 
[dv3, dv4, dt4, O_bt] = bitangentTransfer(O_cper, O_end, 'aa'); % testando tutte e 4 le possibilità conviene aa 

dt5 = TOF(O_end,pi,0); 

DV_1=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4);



% Calcolo il tempo impiegato 

% tempo dalla posizione iniziale fino al punto di cambio di piano
dt1 = TOF(O_start, th_start, th_cp); 

% tempo dalla posizione di cambio di piano alla posizione di cambio di anomalia di pericentro
dt2 = TOF(O_cp, th_cp, th_best(1)); 

% tempo dalla posizione dopo il cambio di anomalia di pericentro fino al pericentro:
dt3 = TOF(O_cper,th_best(1),pi); 

% tempo di manovra della bitangente dt4

% tempo di attesa sull'orbita finale fino al punto d'arrivo
dt5 = TOF(O_end,pi,0); 


DT_1= dt1 + dt2 + dt3 + dt4 + dt5;

dtime = seconds(DT_1);
    dtime.Format = 'hh:mm:ss';

fprintf("Costo della manovra standard: %d \n", DV_1);
fprintf("Tempo impiegato: %s \n", dtime);



Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 1: cambio piano - cambio pericentro - trasferimento bitangente', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,th_cp,dth,'b')
plotOrbit(O_cp,th_cp,th_best(1),dth,'m') % Orbita post cambio piano
plotOrbit(O_cper,th_best(2),pi,dth,'g')
plotOrbit(O_bt,0,pi,dth,'c')
plotOrbit(O_end,pi,0,dth,'r')

% plotOrbit(O_start,0,2*pi,dth, ['--','k'])
% h = findobj(gca, 'Type', 'line', '-not', 'HandleVisibility', 'off');
% h(1).HandleVisibility = 'off';  % Nasconde dalla legenda
% plotOrbit(O_end,0,2*pi,dth,['--','k'])
% h = findobj(gca, 'Type', 'line', '-not', 'HandleVisibility', 'off');
% h(1).HandleVisibility = 'off';  % Nasconde dalla legenda
% xlim([-0.6e5,0.6e5])
% ylim([-0.6e5,0.6e5])
% zlim([-0.6e5,0.6e5])

legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita modificata di piano','Orbita modificata anomalia pericentro','orbita bitangente','Orbita finale')

%%

% sol_1=figure;
% sol_1.Name="Strategia 1: cambio piano - cambio pericentro - trasferimento bitangente";
% scatter3(0,0,0)
% hold on
% scatter3(rr_start(1),rr_start(2),rr_start(3))
% scatter3(rr_end(1),rr_end(2), rr_end(3))
% sist_can(1)
% plotOrbit(O_start,th_start,th_cp,dth,'b') % Orbita di partenza
% plotOrbit(O_cp,th_cp,th_best(1),dth,'r') % Orbita post cambio piano
% plotOrbit(O_cper,th_best(2),pi,dth,'g')
% plotOrbit(O_bt,0,pi,dth,'c')
% plotOrbit(O_end,pi,0,dth,'m')
% xlim([-1e5,1e5])
% ylim([-1e5,1e5])
% zlim([-1e5,1e5])
% legend('attrattore','partenza','arrivo','sitema riferimento','orbita di inizio','orbita modificata di piano','orbita modificata anomalia pericentro','orbita biellittica ausiliaria','orbita finale')

%% Soluzione 2  TEB - CP - CPer 

% Calcolo il costo 

%1) Trasferimento bitangente su orbita ausiliaria 
[dv1, dv2, dt1, O_bt,th0, thf, O_aus] = bitangentTransfer(O_start, O_end, 'aa');
dt2 = TOF(O_start, th_start, pi); 

%2)modifico anomalia pericentro
[dv3, th_cp, O_cp] = changeOrbitalPlane(O_aus, O_end);
dt1 = TOF(O_start, th_start, th_cp);
%3) modifico angolo ppericentro
[dv4, th_i, th_f, th_best, O_cper] = change_pericenter_arg(O_cp, O_end.om, th_cp);
dt2 = TOF(O_cp, th_cp, th_best);
DV_2=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)



Terra_3D(r_terra)
set(gcf, 'Name', 'Strategia 1: cambio piano - cambio pericentro - trasferimento bitangente', 'NumberTitle', 'off');
scatter3(rr_start(1),rr_start(2),rr_start(3))
hold on
scatter3(rr_end(1),rr_end(2), rr_end(3))
plotOrbit(O_start,th_start,pi,dth,'b')
plotOrbit(O_bt,0,pi,dth,'c')
plotOrbit(O_aus,pi,th_cp+2*pi,dth,'c')
plotOrbit(O_cp,th_cp+2*pi,th_best(1),dth,'m')
plotOrbit(O_cper,th_best(2),2*pi,dth,'r')
% plotOrbit(O_end,0,th_end,dth,'y')
% plotOrbit(O_start,0,th_end,dth,'y')
legend('Attrattore','Partenza','Arrivo','Orbita iniziale','Orbita bitangente','Orbita ausiliaria', 'Orbita cambio piano', 'Orbita cambio pericentro', 'Orbita finale')


%% Soluzione   TBT - CP - CPer - TBT modo 1

%1) cambio forma
ra_t = 2 * O_start.a*(1+O_start.e); % n volta l'apocentro di O_start
[dv1, dv2, dv3, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_start,O_end, ra_t);%modifica a,e
%2) cambio piano
[dv4, th_cp, O_cp] = changeOrbitalPlane(O_biel1, O_end);
[dv5, th_i, th_f, th_best, O_cper] = change_pericentre_arg(O_cp, O_end.om, th_cp);
%3) finisco biellittico
O_biel2.i=O_cper.i;
O_biel2.OM=O_cper.OM;
O_biel2.om= O_cper.om;


DV_3=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+ abs(dv5)

fig_3=figure;
fig_3.Name='caso biellittico';
scatter3(0,0,0)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can
plotOrbit(O_start,th_start,2*pi,dth,'b')
plotOrbit(O_biel1,2*pi,th_cp,dth, 'm') 
plotOrbit(O_cp,th_cp,th_best(1),dth, 'c')
plotOrbit(O_cper,th_best(2),pi,dth, 'k')
plotOrbit(O_biel2,pi,th_end, dth, 'g')

plotOrbit(O_end,0,2*pi,dth, 'r')
legend ('attrattore','partenza','arrivo','sitema riferimento','orbita iniziale','prima biell','orbita cambio piano', 'cambio pericentro', 'seconda biell', 'end')

%% Soluzione   TBT - CP - TBT - CPer  modo 2
% questa così non so se è sbagliata o no 
%1) cambio forma
ra_t = 2 * O_start.a*(1+O_start.e); % n volta l'apocentro di O_start
[dv1, dv2, dv3, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_start,O_end, ra_t);%modifica a,e
%2) cambio piano
[dv4, th_cp, O_cp] = changeOrbitalPlane(O_biel1, O_end);
%3) finisco biellittico
O_biel2.i=O_cp.i;
O_biel2.OM=O_cp.OM;
O_biel2.om= O_cp.om;

[dv5, th_i, th_f, th_best, O_cper] = change_pericentre_arg(O_biel2, O_end.om, th_cp);
DV_3=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+ abs(dv5)

 

fig_3=figure;
fig_3.Name='caso biellittico';
scatter3(0,0,0)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can
plotOrbit(O_start,th_start,2*pi,dth,'b')
plotOrbit(O_biel1,2*pi,th_cp,dth, 'm') 
plotOrbit(O_cp,th_cp,pi,dth, 'c')
plotOrbit(O_biel2,pi,th_best(1),dth, 'k') % qui boh perché la faccio concludere in th_best(1) che non è l'apocentro poi la ruoto e il punto di intersazione con end è proprio il punto di arrivo che coincide con il pericentro e quindi faccio qui ultimo impulso
plotOrbit(O_cper,th_best(2),th_end, dth, 'g')
plotOrbit(O_end,0,th_end,dth, 'r')
legend ('attrattore','partenza','arrivo','sitema riferimento','orbita iniziale','prima biell','orbita cambio piano', 'seconda biell', 'cambio pericentro', 'end')

%% Soluzione   TBT - CP - CPer modo 3

%1) cambio forma
ra_t = 2 * O_start.a*(1+O_start.e); % n volta l'apocentro di O_start
[dv1, dv2, dv3, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_start,O_end, ra_t);%modifica a,e
%2) cambio piano
[dv4, th_cp, O_cp] = changeOrbitalPlane(O_biel2, O_end);

[dv5, th_i, th_f, th_best, O_cper] = change_pericentre_arg(O_cp, O_end.om, th_cp);
DV_3=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+ abs(dv5)
 

fig_3=figure;
fig_3.Name='caso biellittico';
scatter3(0,0,0)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can
plotOrbit(O_start,th_start,2*pi,dth,'b')
plotOrbit(O_biel1,2*pi,pi,dth, 'm') 
plotOrbit(O_biel2,pi,th_cp,dth, 'k')
plotOrbit(O_cp,th_cp,th_best(1),dth, 'c')
plotOrbit(O_cper,th_best(2),th_end, dth, 'g')
plotOrbit(O_end,0,th_end,dth, 'r')
legend ('attrattore','partenza','arrivo','sitema riferimento','orbita iniziale','prima biell', 'seconda biell','orbita cambio piano', 'cambio pericentro', 'end')



%% soluzione   CP - TEB - CPer  
% DA CONTROLLARE 
% 1) piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);  
dt1 = TOF(O_start, th_start, th_cp);

[dv3, dv4, dt3, O_bt] = bitangentTransfer(O_cp, O_end, 'aa');

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
plotOrbit(O_start,th_start,th_cp,dth,'b') 
plotOrbit(O_cp,th_cp,pi,dth,'r')
plotOrbit(O_bt,0,th_best(1)-pi,dth,'g')
plotOrbit(O_cper,th_best(2)-pi,pi,dth,'c')
%plotOrbit(O_cper,0,2*pi,dth,'c')
plotOrbit(O_end,0,th_end,dth,'m')
xlim([-1e5,1e5])
ylim([-1e5,1e5])
zlim([-1e5,1e5])
legend('attrattore','partenza','arrivo','sitema riferimento','orbita di inizio','orbita modificata di piano','orbita bitangente','orbita cp','orbita finale')

%% soluzione   CP - TEB - CPer  
% 1) piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);  
dt1 = TOF(O_start, th_start, th_cp);

[dv3, dv4, dt3, O_bt] = bitangentTransfer(O_cp, O_end, 'pa');

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
plotOrbit(O_start,th_start,th_cp,dth,'b') 
plotOrbit(O_cp,th_cp,0,dth,'r')
plotOrbit(O_bt,0,th_best(1),dth,'g')
plotOrbit(O_cper,th_best(2),pi,dth,'c')
plotOrbit(O_end,0,th_end,dth,'m')
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
plotOrbit(O_start,th_start,th_cp,dth,'b') 
plotOrbit(O_cp,th_cp,pi,dth,'r')
plotOrbit(O_bt,pi,th_best(1),dth,'g')
plotOrbit(O_cper,th_best(2),pi+pi,dth,'c')
plotOrbit(O_end,0,th_end,dth,'m')
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
plotOrbit(O_start,th_start,th_cp,dth,'b') 
plotOrbit(O_cp,th_cp,0,dth,'r')
plotOrbit(O_bt,pi + pi,th_best(1),dth,'g')
plotOrbit(O_cper,th_best(2),pi,dth,'c')
plotOrbit(O_end,0,th_end,dth,'m')
xlim([-1e5,1e5])
ylim([-1e5,1e5])
zlim([-1e5,1e5])
legend('attrattore','partenza','arrivo','sitema riferimento','orbita di inizio','orbita modificata di piano','orbita bitangente','orbita cp','orbita finale')

%% soluzione   CP - CPer - TBT
% c'è qualche problema
% 1) piano orbita
[dv1, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end); %la funzione calcola il costo solo in funzione del piano di arrivo, l'orbita che restituisce non è O_end
dt1 = TOF(O_start, th_start, th_cp);
%2)modifico anomalia pericentro
[dv2, thi_cper, thf_cper, th_best, O_cper] = change_pericentre_arg(O_cp, O_end.om-pi, th_cp);
dt2 = TOF(O_cp, th_cp, th_best);
%3)modifico orbita con bitangente 
ra_t = 2 * O_cper.a*(1+O_cper.e);
[dv3, dv4, dv5, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_cper,O_end, ra_t); 

DV_1=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+abs(dv5)

sol_1=figure;
sol_1.Name="soluzione 1: cambio piano - cambio pericentro - trasferimento bitangente";
scatter3(0,0,0)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can
plotOrbit(O_start,th_start,th_cp,dth,'b') % Orbita di partenza
plotOrbit(O_cp,th_cp,th_best(1),dth,'r') % Orbita post cambio piano
plotOrbit(O_cper,th_best(2),0,dth,'g')
plotOrbit(O_biel1,0,pi,dth,'c')
plotOrbit(O_biel2,pi,0,dth,'c')
plotOrbit(O_end,0,th_end,dth,'m')
xlim([-1e5,1e5])
ylim([-1e5,1e5])
zlim([-1e5,1e5])
legend('attrattore','partenza','arrivo','sitema riferimento','orbita di inizio','orbita modificata di piano','orbita modificata anomalia pericentro','orbita biellittica ausiliaria','orbita finale')

%% Soluzione    CP - TBT - CPer

[dv4, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);

ra_t = 2 * O_start.a*(1+O_start.e); % n volta l'apocentro di O_start
[dv1, dv2, dv3, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_cp,O_end, ra_t);%modifica a,e

[dv5, th_i, th_f, th_best, O_cper] = change_pericentre_arg(O_biel2, O_end.om, th_cp);
DV_3=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+ abs(dv5)
 

fig_3=figure;
fig_3.Name='caso biellittico';
scatter3(0,0,0)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can
plotOrbit(O_start,th_start,th_cp,dth,'b')
plotOrbit(O_cp,th_cp,0,dth, 'c')
plotOrbit(O_biel1,0,pi,dth, 'm') 
plotOrbit(O_biel2,pi,th_best(1),dth, 'k')
plotOrbit(O_cper,th_best(2),th_end, dth, 'g')
plotOrbit(O_end,0,th_end,dth, 'r')
legend ('attrattore','partenza','arrivo','sitema riferimento','orbita iniziale','orbita cambio piano','prima biell', 'seconda biell', 'cambio pericentro', 'end')


%% Soluzione    CP - TBT - CPer - TBT

[dv4, th_cp, O_cp] = changeOrbitalPlane(O_start, O_end);

ra_t = 2 * O_start.a*(1+O_start.e); % n volta l'apocentro di O_start
[dv1, dv2, dv3, d_t1, d_t2, O_biel1, O_biel2] = biellipticTransfer(O_cp,O_end, ra_t);%modifica a,e

[dv5, th_i, th_f, th_best, O_cper] = change_pericentre_arg(O_biel1, O_end.om, th_cp);
DV_3=abs(dv1)+abs(dv2)+abs(dv3)+abs(dv4)+ abs(dv5)
 
O_biel2.i=O_cper.i;
O_biel2.OM=O_cper.OM;
O_biel2.om= O_cper.om;

fig_3=figure;
fig_3.Name='caso biellittico';
scatter3(0,0,0)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can
plotOrbit(O_start,th_start,th_cp,dth,'b')
plotOrbit(O_cp,th_cp,0,dth, 'c')
plotOrbit(O_biel1,0,th_best(1),dth, 'm') 
plotOrbit(O_cper,th_best(2),pi, dth, 'g')
plotOrbit(O_biel2,pi,th_end,dth, 'k')
plotOrbit(O_end,0,th_end,dth, 'r')
legend ('attrattore','partenza','arrivo','sitema riferimento','orbita iniziale','orbita cambio piano','prima biell', 'cambio pericentro','seconda biell', 'end')


%%
% valutazione soluzione 3: cambio piano e om nell'orbita di trasferimento 
% della bitangente

% 1) bitangente (i costi di trasferimento saranno uguali)
delta_t1 = TOF(O_start, th_start, pi); % Attesa per la manovra 1
[delta_v1_bt, delta_v2_bt, delta_t4, orbit_bt_temporanea, th0, thf, orbit_forma] = bitangentTransfer(O_start, O_end, 'aa');
% 2) modifico piano
[delta_v2, th_cp, orbit_cp] = changeOrbitalPlane(orbit_bt_temporanea, O_end); % Mi metto nel piano definitivo
delta_t2 = TOF(orbit_bt_temporanea, pi, th_cp); % Attesa per il cambio piano
%3) cambio anomalia del pericentro
[delta_v3, th_i, th_f, th_best, O_cper] = change_pericentre_arg(orbit_cp, O_end.om, th_cp);
delta_t3 = TOF(orbit_cp, th_cp, th_best(1)); % Attesa per cambio orientazione

DELTA_V_3 = abs(delta_v1_bt) + abs(delta_v2_bt) + abs(delta_v2) + abs(delta_v3)
DELTA_T_3 = delta_t1 + delta_t2 + delta_t3 + delta_t4

%plot forma e poi attitudine
fig_2=figure;
fig_2.Name='sol 3: impulso 1 - cambio piano - cambio pericentro - impulso 2';
scatter3(0,0,0)
hold on
scatter3(rr_start(1),rr_start(2),rr_start(3))
scatter3(rr_end(1),rr_end(2), rr_end(3))
sist_can
plotOrbit(O_start, th_start, pi, dth, 'b') %           Start
plotOrbit(orbit_bt_temporanea, 0, th_cp, dth, 'c') %      Trasferimento
plotOrbit(orbit_cp, th_cp, th_best(1), dth, 'g') %        O_cp
plotOrbit(O_cper, th_best(2), pi,dth,'k') %          O_chper
plotOrbit(O_end, pi, th_end,dth,'m') % End
% plotOrbit(orbit_chper,th_best(2),pi,dth,'r')
legend('attrattore','partenza', 'arrivo', 'sistema ref','start','bitangente ausiliaria','orbita di piano','orbita pericentro', 'end')

