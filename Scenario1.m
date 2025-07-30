%% Scenario 1
% Sequenze di manovre standard 
clear all
close all
clc

%parametri figura
dth=0.001;

% Parametro attrattore
mu=398600;
% Dati forniti dell'orbita di partenza
O_start.a = 24400.0;
O_start.e = 0.728300;
O_start.i = 0.104700;
O_start.OM = 1.514000;
O_start.om = 3.107000;
O_start.mu = mu;
th_start = 1.665000;


% Dati forniti dell'orbita d'arrivo
rr=[-12985.280000 3801.011400 8109.619300]';
vv=[-0.948600 -6.134000 1.356000]';

% Trovo dati mancanti dell'orbita d'arrivo 
[O_end,th_end] = car2par(rr,vv,mu);

%valutazione soluzione 1: standard
%1) piano orbita
[delta_v1, om_cp, th_cp, orbit_cp] = changeOrbitalPlane(O_start, O_end); %la funzione calcola il costo solo in funzione del piano di arrivo, l'orbita che restituisce non è O_end
delta_t1 = TOF(O_start, th_start, th_cp);
%2)modifico anomalia pericentro
[delta_v2, th_i_cper, th_f_cper, th_best, orbit_chper] = change_pericentre_arg(orbit_cp, O_end.om-pi, th_cp);
delta_t2 = TOF(orbit_cp, th_cp, th_best);
%3)modifico orbita con bitangente 
[delta_v1_bt, delta_v2_bt, delta_t_bt, orbit_bt] = bitangentTransfer(orbit_chper, O_end, 'aa'); %testanto tutte 4 le possibilità conviene aa con la modifica del pericentro fatta prima 

DELTA_V_1=delta_v1+delta_v2+delta_v1_bt+delta_v2_bt
sol_1=figure;
sol_1.Name="soluzione 1: cambio piano - cambio pericentro - trasferimento bitangente";
scatter3(0,0,0)
hold on
sist_can
plotOrbit(O_start,0,th_cp,dth,'b')
plotOrbit(orbit_cp,th_cp,th_best(1),dth,'r')
plotOrbit(orbit_chper,th_best(2),pi,dth,'g')
plotOrbit(orbit_bt,0,pi,dth,'c')
plotOrbit(O_end,pi,0,dth,'m')
xlim([-1e5,1e5])
ylim([-1e5,1e5])
zlim([-1e5,1e5])
legend('attrattore','sitema riferimento','orbita di inizio','orbita modificata di piano','orbita modificata anomalia pericentro','orbita biellittica ausiliaria','orbita finale')


%valutazione soluzione 2: cambio forma e poi attitudine
%1)modifico orbita con bitangente 

%1) piano orbita
[delta_v1_bt, delta_v2_bt, delta_t, orbit_bt_temporanea, th0, thf, orbit_forma] = bitangentTransfer(O_start, O_end, 'aa');
%2)modifico anomalia pericentro
[delta_v1, om_cp, th_cp, orbit_cp] = changeOrbitalPlane(orbit_forma, O_end);
delta_t1 = TOF(O_start, th_start, th_cp);
%3) modifico angolo ppericentro
[delta_v2, th_i, th_f, th_best, orbit_chper] = change_pericentre_arg(orbit_cp, O_end.om, th_cp);
delta_t2 = TOF(orbit_cp, th_cp, th_best);
DELTA_V_2=delta_v1+delta_v2+delta_v2_bt+delta_v1_bt

%plot forma e poi attitudine
fig_2=figure;
fig_2.Name='caso modifica forma e poi attitudine';
scatter3(0,0,0)
hold on
sist_can
plotOrbit(O_start,0,pi,dth,'b')
plotOrbit(orbit_bt_temporanea,0,pi,dth,'g')
plotOrbit(orbit_forma,pi,th_cp+2*pi,dth,'c')
plotOrbit(orbit_cp,th_cp+2*pi,th_best(1),dth,'m')
plotOrbit(orbit_chper,th_best(2),pi,dth,'r')
legend('attrattore','sistema ref','start','bitangente ausiliaria','orbita di forma','orbita di piano','orbita pericentro')

%valutazione soluzione 3: cambio piano lontano

%1) cambio forma
ra_t = 1 * O_start.a*(1+O_start.e); % n volta l'apocentro di O_start
[delta_v1_be, delta_v2_be, delta_v3_be, delta_t1, delta_t2, orbit_biel1, orbit_biel2] = biellipticTransfer(O_start,O_end, ra_t);%modifica a,e
%2) cambio piano
[delta_v_p, om_cp, th_cp, orbit_cp] = changeOrbitalPlane(orbit_biel1, O_end);
%3) finisco biellittico
orbit_biel2.i=orbit_cp.i;
orbit_biel2.OM=orbit_cp.OM;
orbit_biel2.om= orbit_cp.om;
orbit_biel3=orbit_biel2;
orbit_biel3.e=O_end.e;
orbit_biel3.a=O_end.a;
%4) cambio anomalia pericentro
[delta_v2, th_i, th_f, th_best, orbit_chper] = change_pericentre_arg(orbit_biel2, O_end.om, th_cp);

fig_3=figure;
fig_3.Name='caso biellittico con modifica piano in apocentro';
scatter3(0,0,0)
hold on
plotOrbit(O_start,0,2*pi,dth, 'b')
plotOrbit(orbit_biel1,0,th_cp,dth, 'm') 
plotOrbit(orbit_cp,th_cp,pi,dth, 'c')
plotOrbit(orbit_biel2,pi,0,dth, 'k')
plotOrbit(orbit_biel3,0,th_best(1),dth, 'g')

plotOrbit(O_end,0,2*pi,dth, 'r')
legend ('terra','start','prima biell','transfer al cambio piano','seconda biell','transfer al pericentro','end')

%COMMENTO: da valutare se ha senso modificare due volte il pericentro
%facendo si che il cambio piano avvenga da biellittica 1 direttamente a
%biellittica 2







% Supponendo che delta_v2 rimane uguale dopo al cambio piano
DELTA_V_sol3 = delta_v1_be + delta_v_p + delta_v2_be + delta_v3_be + delta_v2
% con ra_t 10*ra_start, 6.6643
% con ra_t 20*ra_start, 6.8317 :(
% con ra_t 5*ra_start, 6.3390 :)
% con ra_t 4*ra_start, 6.1810 :|
% con ra_t 1*ra_start, 4.6866 XOXO tragedia

% prova plot soluzione prima forma poi attitude con bitangente aa
% th0=0;
% thf=2*pi;
% dth=pi/100;
% figure
% scatter3(0,0,0)
% hold on
% plotOrbit(O_start,th0,thf,dth, 'b')
% plotOrbit(orbit_forma,th0,thf,dth, 'y')
% plotOrbit(orbit_cp,th0,thf,dth, 'p')
% plotOrbit(O_end,th0,thf,dth, 'r')
% legend('attrattore','start','forma','piano','finale')
%%







% Plotto le due orbite
th0=0;
thf=2*pi;
dth=pi/100;

plotOrbit(O_start,th0,thf,dth, 'b') % Partenza
hold on
plotOrbit(O_end,th0,thf,dth, 'k') % Arrivo

% Plotto il pianeta Terra
[X,Y,Z]=sphere(1000);
k=0.5*10^4;
surf(k*X,k*Y,k*Z,"LineStyle","none","FaceColor","texturemap",CData=flip((imread('earth.jpg'))));
xlim([-0.6*10^5;0.6*10^5])
ylim([-0.6*10^5;0.6*10^5])
zlim([-0.6*10^5;0.6*10^5])

% Modifico il piano dell'orbita iniziale 
[delta_v1, om_cp, th_cp, orbit_cp] = changeOrbitalPlane(O_start, O_end);
plotOrbit(orbit_cp,th0,thf,dth,['--', 'r'])
delta_t1 = TOF(O_start, th_start, th_cp);

% Modifico anomalia del pericentro
[delta_v2, th_i, th_f, th_best, orbit_chper] = change_pericentre_arg(orbit_cp, O_end.om, th_cp);
plotOrbit(orbit_chper,th0,thf,dth,['--', 'g'])
delta_t2 = TOF(orbit_cp, th_cp, th_best);

delta_t3 = TOF(orbit_cp, th_cp, th_best);
% Effettuo il trasferimento di orbita 
% [delta_v1_bt, delta_v2_bt, delta_t_bt, orbit_bt] = bitangentTransfer(orbit_chper, O_end, 'ap');
% plotOrbit(orbit_bt,th0,pi,dth,['-', 'c'])
% non penso abbiano molto senso i trasferimenti aa o pp per come sono messe
% le orbite 




% il raggio di apocentro dell'orbita iniziale è 42171, ra_t è arbitrario
ra_t = 45000.0;
[delta_v11, delta_v22, delta_v33, delta_t11, delta_t22, orbit_biel1, orbit_biel2] = biellipticTransfer(orbit_chper,O_end, ra_t);
plotOrbit(orbit_biel1,th0,pi,dth,['-', 'c'])
hold on
plotOrbit(orbit_biel2,pi,2*pi,dth,['-', 'r'])


