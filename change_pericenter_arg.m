function [delta_v, th_i, th_f, th_best, orbit_chper] = change_pericenter_arg(orbit, om_f, th_0)
% La funzione prende in ingresso il semiasse maggiore, l'eccentricità,
% l'anomalia del pericentro iniziale e finale e l'anomalia vera di partenza.
% Restituisce il costo della manovra di cambio dell'anomalia del pericentro, 
% le anomalie vere alle quali può verificarsi la manovra, l'anomalia vera 
% migliore per effettuare la manovra e i parametri completi dell'orbita d'arrivo.


% Input
% orbit struttura contenente i parametri fondamentali dell'orbita:
%  orbit.a   valore del semiasse maggiore
%  orbit.e   modulo del vettore eccentricità
%  orbit.om  anomalia del pericentro iniziale
%  orbit.mu  parametro gravitazionale dell'attrattore in questione
%  om_f      anomalia del pericentro finale
%  th_0      anomalia vera di partenza

% Output
%  delta_v    costo manovra  
%  th_i       vettore [2x1] contenente le anomalie vere di manovra in
%             riferimento all'orbita iniziale 
%  th_f       vettore [2x1] contenente le anomalie vere di manovra in
%             riferimento all'orbita finale 
%  th_best    anomalia vera migliore dove effettuare la manovra in
%             riferimento all'orbita iniziale e all'orbita finale
% orbit_chper è una struttura che contiene i parammetri caratterizzanti
% dell'orbita finale dopo il cambio di anomalia del pericentro

%  orbit_chper.a   valore del semiasse maggiore iniziale e finale 
%  orbit_chper.e   modulo del vettore eccentricità iniziale e finale 
%  orbit_chper.i   angolo d'inclinazione dell'orbita iniziale e finale
%  orbit_chper.OM  ascensione retta del nodo ascendete (RAAN) iniziale e finale
%  orbit_chper.om  anomalia del pericentro finale 

% Parametri oribita di partenza
a = orbit.a;
e = orbit.e;
om_i = orbit.om;
mu = orbit.mu;
D_om = om_f - om_i; % variazione anomalia del pericentro

% Calcolo orbita finale
orbit_chper = orbit;
orbit_chper.om = om_f;

% theta di manovra in riferimento iniziale
th_1i = D_om/2;
th_2i = D_om/2+pi;
th_i = [th_1i;th_2i];

% theta di manovra in riferimento d'arrivo
th_1f = 2*pi-D_om/2;
th_2f = pi-D_om/2;
th_f = [th_1f;th_2f];


p=a*(1-e^2); % semilato retto
delta_v = 2*sqrt(mu/p)*e*sin(abs(D_om/2)); % costo della manovra 

if D_om <= 0
    disp('Delta omega nel trasferimento di pericentro è negativo')
end
th_best=zeros(2,1);
% Cerco il theta migliore; th_0 è nel primo sistema di riferimento, i punti
% d'intersezione sono quindi D_om/2 e D_om/2+pi. Per confrontarli, ho vari
% casi:
if th_0 == th_i(1) 
    th_best(1) = th_i(1);
    th_best(2) = th_f(1);
elseif th_0 == th_i(2)
        th_best(1) = th_i(2);
        th_best(2) = th_f(2);
elseif th_0 < 0 
    if D_om < 0
        if th_0 > th_i(1) % Quindi sono tra th1 e th2, vado in th2
            th_best(1) = th_i(2); % Il th2 del sdr 1
            th_best(2) = th_f(2); % Il th2 del sdr 2
        else % Sono tra th2 e th1, vado in th1
            th_best(1) = th_i(1);
            th_best(2) = th_f(1);
        end
    else % D_om > 0 e th_0 < 0
        th_0 = th_0 + 2*pi; % Quindi th0, th1, th2 tutti compresi tra 0 e 2pi
        if th_0 > th_i(1) && th_0 < th_i(2) % Sono tra th1 e th2, vado in th2
            th_best(1) = th_i(2); % Il th2 del sdr 1
            th_best(2) = th_f(2); % Il th2 del sdr 2
        else % Sono tra th2 e th1, vado in th1
            th_best(1) = th_i(1);
            th_best(2) = th_f(1);
        end    
    end
else % th_0 > 0 % Li ho tutti e tre compresi tra 0 e 2pi
    if th_0 > th_i(1) && th_0 < th_i(2) % Vado in th2
        th_best(1) = th_i(2);
        th_best(2) = th_f(2);
    else % Vado in th1
        th_best(1) = th_i(1);
        th_best(2) = th_f(1);
    end
end 


end
