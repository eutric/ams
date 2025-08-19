function [rr,vv, matrici] = par2car(orbit, th)
% La funzione prende in ingresso i parametri caratterizzanti nel sistema di
% riferimento perifocale (Kepleriano) e restituisce i vettori posizione e 
% velocità di un punto dell'orbita in coordinate cartesiane (sistema di 
% riferimento Geocentrico inerziale)

% Input
% orbit struttura contenente tutti i parametri fondamentali dell'orbita:
%  orbit.a   valore del semiasse maggiore
%  orbit.e   modulo del vettore eccentricità
%  orbit.i   angolo d'inclinazione dell'orbita
%  orbit.OM  ascensione retta del nodo ascendente (RAAN)
%  orbit.om  anomalia del pericentro
%  orbit.mu  Parametro gravitazionale dell'attrattore in questione
% th anomalia vera

% Output
% rr  Vettore [3x1] contenente le coordinate assolute della posizione 
% vv  Vettore [3x1] contenente le coordinate assolute della velocità

a=orbit.a;
e=orbit.e;
i=orbit.i;
om=orbit.om;
OM=orbit.OM;
mu=orbit.mu;

% Rotazione di OM intorno al versore k
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

p=a*(1-e^2); % semilato retto
r=p/(1+e*cos(th));  % modulo della posizione
rr2=r*[cos(th); sin(th);0];%  % vettore posizione nel perifocale (2 xke in 2D)
vv2=sqrt(mu/p)*[-sin(th);e+cos(th);0]; % vettore velocità nel perifocale

T=R_om*R_i*R_OM; % matrice di trasformazione ECI -> PF
% T' PF -> ECI
rr=T'*rr2; % vettore posizione nel cartesiano 
vv=T'*vv2; % vettore velocità nel cartesiano 

matrici.R_om = R_om;
matrici.R_i = R_i;
matrici.R_OM = R_OM;
matrici.T = T;
end

