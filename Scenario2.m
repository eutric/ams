clear
close all
clc

mu = 1.32712440042e20 * 0.001^3; % km^3/s^2
dth = 0.001;

% Terra
O_start.a = 1.4946e8;  % [km]
O_start.e = 0.016;     % [ ]
O_start.i = 9.1920e-5; % [rad]
O_start.OM = 2.7847;   % [rad]
O_start.om = 5.2643;   % [rad]
O_start.mu = mu;       % [km^3/s^2]

% Asteroide 163899 (2003 SD220)
O_end.a = 0.827903 *1.496e8; % [km]
O_end.e = 0.209487;          % [ ]
O_end.i = 8.55 *pi/180;      % [rad]
O_end.OM = 273.63 *pi/180;   % [rad]
O_end.om = 327.03 *pi/180;   % [rad]
O_end.mu = mu;               % [km^3/s^2]

figure
scatter3 (0,0,0, 'red', LineWidth=2)
hold on
plotOrbit (O_start, 0, 2*pi, dth, 'b')
plotOrbit (O_end, 0, 2*pi, dth, 'k')
%prova caso esempio scelgo theta1i=theta2f=0
th1i=0;
th2f=0;
[rri,vvi] = par2car(O_start, th1i);
[rrf,vvf] = par2car(O_end, th2f);
scatter3(rri(1),rri(2),rri(3))
scatter3(rrf(1),rrf(2),rrf(3))

%identifico orbita
hh=cross(rri,rrf)./(norm(cross(rri,rrf)));
i=acos(dot(hh,[0,0,1]));
NN=cross([0,0,1],hh)./norm(cross([0,0,1],hh));

if dot(NN,[0,1,0])>=0
    OM=acos(dot(NN,[1,0,0]));
elseif dot(NN,[0,1,0])<0
    OM=2*pi-acos(dot(NN,[1,0,0]));
end
R_OM=[cos(OM) sin(OM) 0;
    -sin(OM) cos(OM) 0;
    0 0 1];
%scelgo 
om=0;
% Rotazione di i intorno al versore i'
R_i=[1 0 0;
    0 cos(i) sin(i);
    0 -sin(i) cos(i)];

% Rotaziozione di om intorno al versore k''
R_om=[cos(om) sin(om) 0;
    -sin(om) cos(om) 0;
    0 0 1];
T_e_pf=R_om*R_i*R_OM;
rri_t=T_e_pf*rri;
rrf_t=T_e_pf*rrf;
th1t=atan(rri(2)/rri(1));
th2t=atan(rrf(2)/rrf(1));
if th1t<=0 || th2t<=0
    th1t=th1t+2*pi;
    th2t=th2t+2*pi;
end
%calcolo valori scalari e parametri di forma
r1=norm(rri_t);
r2=norm(rrf_t);
et=(r2-r1)/(r1*cos(th1t)-r2*cos(th2t));
at=r1*(1+et*cos(th1t))/(1-et^2);
O_t.a=at;
O_t.e=et;
O_t.i=i;
O_t.OM=OM;
O_t.om=om;
O_t.mu=O_start.mu;
[rr1t,vv1t]=par2car(O_t,th1t);
[rr2t,vv2t]=par2car(O_t,th2t);
cost=abs(norm(vv1t-vvi))+abs(norm(vv2t-vvf));
plotOrbit(O_t,0,2*pi,0.01)
xlim([-2e8,2e8])
ylim([-2e8,2e8])
zlim([-2e8,2e8])

%% test funzione
clear 
close all
clc

mu = 1.32712440042e20 * 0.001^3; % km^3/s^2
dth = 0.001;

% Terra
O_start.a = 1.4946e8;  % [km]
O_start.e = 0.016;     % [ ]
O_start.i = 9.1920e-5; % [rad]
O_start.OM = 2.7847;   % [rad]
O_start.om = 5.2643;   % [rad]
O_start.mu = mu;       % [km^3/s^2]

% Asteroide 163899 (2003 SD220)
O_end.a = 0.827903 *1.496e8; % [km]
O_end.e = 0.209487;          % [ ]
O_end.i = 8.55 *pi/180;      % [rad]
O_end.OM = 273.63 *pi/180;   % [rad]
O_end.om = 327.03 *pi/180;   % [rad]
O_end.mu = mu;          

n = 100; % Dimensioni griglia
a = 0;
b = 2*pi;
OO_t = []; % Vettore di oggetti orbita, non si può preallocare
kk = 0; % Contatore
costmin=1e16;

figure
scatter3 (0,0,0, 'red', LineWidth=2)
hold on
plotOrbit (O_start, 0, 2*pi, dth, 'b')
plotOrbit (O_end, 0, 2*pi, dth, 'k')
tic
for th1i=linspace(a, b, n)
    for th2f=linspace(a, b, n)
        for om=linspace(a, b,n)
            [O_t,th_t,cost] = O_tfun(O_start,O_end,th1i,th2f,om);
            if O_t.e<=1 && O_t.e>=0
                OO_t=[OO_t,O_t];
                % plotOrbit (O_t, 0, 2*pi, dth, 'y--')
                if cost<costmin
                    costmin=cost;
                    O_best=O_t;
                    th_t_best=th_t;
                    th_best_12 = [th1i, th2f];
                end
            else 
                kk = kk+1;
            end

        end
    end
end
% 
[punto_1, ~] = par2car(O_start, th_best_12(1));
[punto1_t, ~] = par2car(O_best, th_t_best(1));
[punto2_t, ~] = par2car(O_best, th_t_best(2));
[punto_2, ~] = par2car(O_end, th_best_12(2));

scatter3 (punto1_t(1), punto1_t(2), punto1_t(3), 'red')
scatter3 (punto2_t(1), punto2_t(2), punto2_t(3), 'red')
scatter3 (punto_1(1), punto_1(2), punto_1(3), 'red')
scatter3 (punto_2(1), punto_2(2), punto_2(3), 'red')
plotOrbit(O_best, th_t_best(1), th_t_best(2), dth, 'g') % Trasferimento vincente

th_best_12
costmin
% th1_t = 3.6368, th2_t = 3.5000, dV = 6.3194