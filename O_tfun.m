function [O_t] = O_tfun(O_start,O_end,th1i,th2f,om)
%proviamo a fare una funzione che nserendo delle limitazioni in input
%restituisce l'orbita di trasferimento in output
[rri,vvi] = par2car(O_start, th1i);
[rrf,vvf] = par2car(O_end, th2f);
%identifico orbita
hh=cross(rri,rrf)./(norm(cross(rri,rrf)));
i=acos(dot(hh,[0,0,1]));
NN=cross([0,0,1],hh)./norm(cross([0,0,1],hh));
if dot(NN,[0,1,0])>=0
    OM=acos(dot(NN,[1,0,0]));
else
    OM=2*pi-acos(dot(NN,[1,0,0]));
end

R_OM=[cos(OM) sin(OM) 0;
    -sin(OM) cos(OM) 0;
    0 0 1];

% Rotazione di i intorno al versore i'
R_i=[1 0 0;
    0 cos(i) sin(i);
    0 -sin(i) cos(i)];

% Rotaziozione di om intorno al versore k''
R_om=[cos(om) sin(om) 0;
    -sin(om) cos(om) 0;
    0 0 1];
T_e_pf=R_om*R_i*R_OM;

rri_t = T_e_pf*rri; % r(theta) nell'orbita di trasferimento
rrf_t = T_e_pf*rrf; % r(theta) nell'orbita di trasferimento
r1 = norm(rri);
r2 = norm(rrf);

cos_th1 = rri_t(1)/r1;
cos_th2 = rrf_t(1)/r2;

if rri_t(2) >= 0
    th1_t = acos (cos_th1);
else
    th1_t = 2*pi - acos (cos_th1);
end

if rrf_t(2)  >= 0
    th2_t = acos(cos_th2);
else
    th2_t = 2*pi - acos(cos_th2);
end

% th1_t = atan(rri_t(2)/rri_t(1)); % Possono esserci problemi con la periodicità?
% th2_t = atan(rrf_t(2)/rrf_t(1));


%calcolo valori scalari e parametri di forma

e_t= (r2-r1)/(r1*cos_th1 - r2*cos_th2);
a_t= r1*(1+e_t*cos_th1)/(1-e_t^2);
O_t.a=a_t;
O_t.e = e_t;
O_t.i=i;
O_t.OM=OM;
O_t.om=om;
O_t.mu=O_start.mu;

if e_t >= 0 && e_t < 1
    [~,vv1t]=par2car(O_t,th1_t);
    [~,vv2t]=par2car(O_t,th2_t);
    O_t.cost1=norm(vv1t-vvi);
    O_t.cost2=norm(vv2t-vvf);
    O_t.cost= norm(vv1t-vvi) + norm(vv2t-vvf);
    O_t.tempo = TOF (O_t, th1_t, th2_t);
else
    O_t.cost = 1e15;
    O_t.tempo = 1e20;
end

O_t.th_t=[th1_t,th2_t]; % th1 e th2 nel perifocale dell'orbita di trasferimento
end





