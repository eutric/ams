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
% ra e rp sono invertitii
ra_i=O_start.a*(1-O_start.e);
rp_i=O_start.a*(1+O_start.e);

%dati orbita arrivo in forma parametrica
rr=[-12985.280000 3801.011400 8109.619300]';
vv=[-0.948600 -6.134000 1.356000]';
[O_end,th] = car2par(rr,vv,mu);
ra_f=O_end.a*(1-O_end.e);
rp_f=O_end.a*(1+O_end.e);
rp_f/rp_i
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
    [delta_v1_bt, delta_v2_bt, delta_t_bt] = bitangentTransfer(O_pp, O_end, type{:});
    cost=delta_v1_bt+delta_v2_bt;
    if cost<min || min==0
        min=cost;
        type_min=type;
    end
end
ra_t_vect=linspace(ra_i,100*O_end.a,1000);
min_bet=0;
ii=1;
cost_bet=zeros(length(ra_t_vect),1);
for ra_ii=ra_t_vect
  
    [delta_v1_bet, delta_v2_bet, delta_v3_bet, delta_t1_bet, delta_t2_bet] = biellipticTransfer(O_pp,O_end, ra_ii);
    delta_v=delta_v2_bet+delta_v3_bet+delta_v1_bet;
    if delta_v<min_bet || min_bet==0
        min_bet=delta_v;
        ra_best=ra_ii;
    end
    cost_bet(ii)=delta_v;  
    ii=ii+1;

        
end
figure
semilogy(ra_t_vect,cost_bet')
hold on
grid on
semilogy(ra_t_vect,min*ones(length(ra_t_vect),1))







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
[rr_i,vv_i] = par2car(O_start, th_start);
% Orbita d'Arrivo
rr=[-12985.280000 3801.011400 8109.619300]';
vv=[-0.948600 -6.134000 1.356000]';
% 
[O_end,th_end] = car2par(rr,vv,mu);

%cambio piano
dth=pi/100;
[delta_v1,om_f, thetacp, o_plane] = changeOrbitalPlane(O_start, O_end);
delta_t1 = TOF(O_start, th_start, thetacp);
figure

% Plot punto di partenza
scatter3(rr_i(1),rr_i(2),rr_i(3),'*',LineWidth=2)
hold on
% Plot punto finale
scatter3(rr(1),rr(2),rr(3),'*',LineWidth=2)
% Plot corpo attrattore
scatter3(0,0,0,LineWidth=5)

%plot prima orbita fino al punto in cui cambio piano
plotOrbit(O_start,th_start,thetacp,dth,'--')

%cambio pericentro
[delta_v2,th_1,th_2,th_best,o_pericentre] = change_pericentre_arg(o_plane,O_end.om,thetacp);
delta_t2 = TOF(o_plane, thetacp, th_best(1));
%plotto orbita a piano modificato fino alla anomalia vera ottimale per
%effettuare il cambio di om
plotOrbit(o_plane,thetacp,th_best(1),dth,'-')

%calcolo teta finale dell orbita di pericentro cambiato, essendo di tipo
%pa, cambio orbita nel pericentro ovvero per teta = 0

[delta_v3_1, delta_v3_2, delta_t4, orbit_bt,th0,thf] = bitangentTransfer(o_pericentre, O_end, ['pa']);
delta_t3 = TOF(o_pericentre, th_best(2), th0(1));
delta_t5 = TOF(O_end,thf(2),th_end);

plotOrbit(o_pericentre, th_best(2), th0(1), dth,'-')
plotOrbit(orbit_bt, th0(3), thf(3), dth, '-') % L'orbita intermedia, la dovrei sempre 
% percorrere da 0 a pi !!! Controllare che vale anche se rimpicciolisco 
% l'orbita
plotOrbit(O_end,thf(2),th_end,dth,'-')
legend('pt partenza','pt arrivo','Corpo Attrattore','orbita di partenza fino a prima manovra','orbita piano cambiato','orbita om cambiato','orbita intermedia','orbita finale fino a teta')

costo=delta_v1+delta_v2+delta_v3_1+delta_v3_2
tempo=delta_t5+delta_t3+delta_t4+delta_t2+delta_t1
tempo_ore = tempo / 3600

% Rapporto da verificare per trasferimento bitangente(Honmann Generalizzato)/biellittico
% tra due ellissi è apocentro_finale/pericentro_iniziale
r_i = o_pericentre.a*(1-o_pericentre.e);
r_f = O_end.a*(1+o_pericentre.e);
rapporto = r_f/r_i


%% STRATEGIA STANDARD 
% Esempio da apocentro a apocentro  
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
[rr_i,vv_i] = par2car(O_start, th_start);
% Orbita d'Arrivo
rr=[-12985.280000 3801.011400 8109.619300]';
vv=[-0.948600 -6.134000 1.356000]';
% 
[O_end,th_end] = car2par(rr,vv,mu);

%cambio piano
dth=pi/100;
[delta_v1,om_f, thetacp, o_plane] = changeOrbitalPlane(O_start, O_end);
delta_t1 = TOF(O_start, th_start, thetacp);
figure

% Plot punto di partenza
scatter3(rr_i(1),rr_i(2),rr_i(3),'*',LineWidth=2)
hold on
% Plot punto finale
scatter3(rr(1),rr(2),rr(3),'*',LineWidth=2)
% Plot corpo attrattore
scatter3(0,0,0,LineWidth=5)

%plot prima orbita fino al punto in cui cambio piano
plotOrbit(O_start,th_start,thetacp,dth,'--')

%cambio pericentro
[delta_v2, th_1, th_2, th_best, o_pericentre] = change_per2(o_plane, O_end.om, thetacp);
delta_t2 = TOF(o_plane, thetacp, th_best(1));
%plotto orbita a piano modificato fino alla anomalia vera ottimale per
%effettuare il cambio di om
plotOrbit(o_plane,thetacp,th_best(1) - 2*pi,dth,'-')

%calcolo teta finale dell orbita di pericentro cambiato, essendo di tipo
%pa, cambio orbita nel pericentro ovvero per teta = 0

[delta_v3_1, delta_v3_2, delta_t4, orbit_bt,th0,thf] = bitangentTransfer(o_pericentre, O_end, ['aa']);
delta_t3 = TOF(o_pericentre, th_best(2), th0(1));
delta_t5 = TOF(O_end,thf(2),th_end);

plotOrbit(o_pericentre, th_best(2), th0(1), dth,'-')
plotOrbit(orbit_bt, th0(3), thf(3), dth, '-') % L'orbita intermedia, la dovrei sempre 
% percorrere da 0 a pi !!! Controllare che vale anche se rimpicciolisco 
% l'orbita
plotOrbit(O_end,thf(2),th_end,dth,'-')
legend('pt partenza','pt arrivo','Corpo Attrattore','orbita di partenza fino a prima manovra','orbita piano cambiato','orbita om cambiato','orbita intermedia','orbita finale fino a teta')

costo=delta_v1+delta_v2+delta_v3_1+delta_v3_2
tempo=delta_t5+delta_t3+delta_t4+delta_t2+delta_t1
tempo_ore = tempo / 3600

% Rapporto da verificare per trasferimento bitangente(Honmann Generalizzato)/biellittico
% tra due ellissi è apocentro_finale/pericentro_iniziale
r_i = o_pericentre.a*(1-o_pericentre.e);
r_f = O_end.a*(1+o_pericentre.e);
rapporto = r_f/r_i


%% Appunti scenario 2
% questo giro è un trasferimento diretto
% Si ignora che i pianeti si muovono durante il traserimento (quindi la
% vera posizione dei corpi attrattori)
