function [] = plotOrbit_plane(orbit,th0,thf,dth,linestyle)
% La funzione prende in ingresso i parametri di un'orbita, l'anomalia vera
% iniziale e finale, la step size e una stringa o un vettore di stringhe e
% restituisce il plot bidimensionale dell'orbita.

% Input
% orbit struttura contenente tutti i parametri fondamentali dell'orbita:
%  orbit.a   valore del semiasse maggiore
%  orbit.e   modulo del vettore eccentricità
%  orbit.i   angolo d'inclinazione dell'orbita 
%  orbit.OM  ascensione retta del nodo ascendete (RAAN) 
%  orbit.om  anomalia del pericentro
%  orbit.mu  parametro gravitazionale dell'attrattore in questione
%  th0       anomalia vera iniziale
%  thf       anomalia vera finale
%  dth       step size

% Output
%  Grafico tridimensionale dell'orbita


% Parametri dell'orbita
a = orbit.a;
e = orbit.e;
i = orbit.i;
om = orbit.om;
OM = orbit.OM;
% mu = orbit.mu; % non lo abbiamo usato ma penso serva per il VETTORE VELOCITA'
% che per ora non abbiamo inserito
if thf<=th0
    thf=thf+2*pi;
end

theta = th0:dth:thf; % angolo discretizzato 
p = a*(1-e^2);  % semilato retto
r = @(t)p./(1+e.*cos(t)); % funzione della posizione nell'orbita
rv = r(theta);
rx = rv.*cos(theta);
ry = rv.*sin(theta);
r2 = [rx; ry; zeros(length(rx),1)'];

rp = [r(0)*cos(0); r(0)*sin(0);0];
ra = [r(pi)*cos(pi); r(pi)*sin(pi);0];

% Matrici di rotazione
R_OM = [cos(OM) sin(OM) 0;
    -sin(OM) cos(OM) 0;
    0 0 1];
R_i = [1 0 0;
    0 cos(i) sin(i);
    0 -sin(i) cos(i)];
R_om = [cos(om) sin(om) 0;
    -sin(om) cos(om) 0;
    0 0 1];

T = R_om*R_i*R_OM; % matrice di trasformazione 
rr = T' * r2; 

rp_car = T' *rp;
ra_car = T' *ra;

if nargin==5 
    plot(rr(1,:),rr(2,:),linestyle,LineWidth=1.4)
else
    plot(rr(1,:),rr(2,:))
end
% hold on
% plot ([rp_car(1), ra_car(1)], [rp_car(2), ra_car(2)], linestyle)
grid on
end

