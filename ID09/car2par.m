function [orbit,th] = car2par(rr,vv,mu)
% La funzione prende in ingresso i vettori posizione e velocità di un punto
% dell'orbita in coordinate cartesiane (sistema di riferimento Geocentrico
% inerziale) e restituisce i parametri caratterizzanti nel sistema di
% riferimento perifocale

% Input
% rr  Vettore [3x1] contenente le coordinate assolute della posizione 
% vv  Vettore [3x1] contenente le coordinate assolute della velocità
% mu  Parametro gravitazionale dell'attrattore in questione

% Output
% orbit struttura contenente tutti i parametri fondamentali dell'orbita:
%  orbit.a   valore del semiasse maggiore
%  orbit.e   modulo del vettore eccentricità
%  orbit.i   angolo d'inclinazione dell'orbita
%  orbit.OM  ascensione retta del nodo ascendete (RAAN)
%  orbit.om  anomalia del pericentro
%  orbit.mu  Parametro gravitazionale dell'attrattore in questione
% th anomalia vera

kk=[   % Versore k, del sistema di riferimento cartesiano
    0;
    0;
    1
];

r=norm(rr); % Modulo del vettore posizione
v=norm(vv); % Modulo del vettore velocità
orbit.a=(2/r-v.^2/mu).^-1; % calcolo semiasse maggiore

hh=cross(rr,vv); % Vettore momento della quantità di moto per unità di massa
h=norm(hh);      % Modulo vettore momento della quantità di moto

ee=cross(vv,hh)./mu-rr./r; % Vettore eccentricità
orbit.e=norm(ee); % Modulo vettore eccentricità

orbit.i=acos(hh(3)./h); % Angolo di inclinazione dell'orbita
N=cross(kk,hh)./norm(cross(kk,hh)); % Versore nodo ascendente
if N(2)>=0
    orbit.OM=acos(N(1)); 
else
    orbit.OM=2*pi-acos(N(1)); % RAAN
end

if ee(3)>=0
    orbit.om=acos(dot(N,ee)/orbit.e);
else
    orbit.om=2*pi-acos(dot(N,ee)/orbit.e); % Anomalia del pericentro
end
vr=dot(vv,rr)/r; % Prodotto scalare, tra posizione e velocità
if vr>=0
    th=acos(dot(rr,ee)/(r*orbit.e));
else
    th=2*pi-acos(dot(rr,ee)/(r*orbit.e));
end

orbit.mu=mu;
end

