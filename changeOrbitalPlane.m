function [delta_v, om_f, theta, orbit_cp] = changeOrbitalPlane(orbit_i, orbit_f)
% La funzione prende in ingresso i parametri orbitali iniziali, l'inclinazione e 
% l'ascensione retta del nodo ascendente dell'orbita finale e restituisce
% il costo della manovra, l'anomalia del pericentro finale, l'anomalia vera
% alla manovra e i parametri completi dell'orbita d'arrivo.

% Input
% orbit_i e orbit_f sono strutture contenenti tutti o alcuni parametri
% fondamentali dell'orbita iniziale o finale dopo il cambio di piano:
%  orbit_i.a   valore del semiasse maggiore iniziale e finale 
%  orbit_i.e   modulo del vettore eccentricità iniziale e finale 
%  orbit_i.i   angolo d'inclinazione dell'orbita iniziale
%  orbit_i.OM  ascensione retta del nodo ascendete (RAAN) iniziale
%  orbit_i.om  anomalia del pericentro iniziale
%  orbit_i.mu  Parametro gravitazionale dell'attrattore in questione
%  orbit_f.OM  ascensione retta del nodo ascendete (RAAN) finale 
%  orbit_f.i   angolo d'inclinazione dell'orbita finale

% Output
%  delta_v      costo della manovra 
%  om_f         anomalia del pericentro finale
%  theta        anomalia vera alla manovra
% orbit_cp è una struttura che contiene i parammetri caratterizzanti
% dell'orbita finale dopo il cambio di piano
%  orbit_cp.a   valore del semiasse maggiore iniziale e finale 
%  orbit_cp.e   modulo del vettore eccentricità iniziale e finale 
%  orbit_cp.i   angolo d'inclinazione dell'orbita finale
%  orbit_cp.OM  ascensione retta del nodo ascendete (RAAN) finale

% Parametri iniziali
a = orbit_i.a;
e = orbit_i.e;
i_i = orbit_i.i;
OM_i = orbit_i.OM;
om_i = orbit_i.om;
mu = orbit_i.mu;

% Parametri finali
i_f=orbit_f.i;
OM_f=orbit_f.OM;

% Creo l'orbita d'arrivo
orbit_cp=orbit_i;
orbit_cp.OM=OM_f;
orbit_cp.i=i_f;


% Risoluzione del triangolo sferico con le varie casistiche 
delta_OM = OM_f - OM_i; % variazione di RAAN
delta_i = i_f - i_i; % variazione dell'inclinazione 
p = a*(1-e^2); % semilato retto

% aggiungi caso in cui delta_OM=0 che il triagolo collassa in un punto e
% delta_i = 0

if delta_OM>0 && delta_i>0
    alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(abs(delta_OM)));
    cosu_i = (-cos(i_f) + cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
    cosu_f = (cos(i_i) - cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
    sinu_i = (sin(abs(delta_OM))*sin(i_f)) / sin(alpha);
    sinu_f = (sin(abs(delta_OM))*sin(i_i)) / sin(alpha);
    u_i = atan2(sinu_i, cosu_i);
    u_f = atan2(sinu_f, cosu_f);

    theta = u_i - om_i;
    om_f = u_f - theta;

    if cos(theta)>0 
        v_theta = sqrt(mu/p)*(1+e*cos(theta + pi));
    else 
        v_theta = sqrt(mu/p)*(1+e*cos(theta));
    end

    delta_v = 2*v_theta*sin(alpha/2);


elseif delta_OM>0 && delta_i<0
    alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(abs(delta_OM)));
    cosu_i = (cos(i_f) - cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
    cosu_f = (-cos(i_i) + cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
    sinu_i = (sin(abs(delta_OM))*sin(i_f)) / sin(alpha);
    sinu_f = (sin(abs(delta_OM))*sin(i_i)) / sin(alpha);
    u_i = atan2(sinu_i, cosu_i);
    u_f = atan2(sinu_f, cosu_f);

    theta = 2*pi - u_i - om_i;
    om_f = 2*pi - u_f - theta;

    if cos(theta)>0 
        v_theta = sqrt(mu/p)*(1+e*cos(theta + pi));
    else 
        v_theta = sqrt(mu/p)*(1+e*cos(theta));
    end

    delta_v = 2*v_theta*sin(alpha/2);

elseif delta_OM<0 && delta_i>0
    alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(abs(delta_OM)));
    cosu_i = (-cos(i_f) + cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
    cosu_f = (cos(i_i) - cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
    sinu_i = (sin(abs(delta_OM))*sin(i_f)) / sin(alpha);
    sinu_f = (sin(abs(delta_OM))*sin(i_i)) / sin(alpha);
    u_i = atan2(sinu_i, cosu_i);
    u_f = atan2(sinu_f, cosu_f);

    theta = 2*pi - u_i - om_i;
    om_f = 2*pi - u_f - theta;

     if cos(theta)>0 
        v_theta = sqrt(mu/p)*(1+e*cos(theta + pi));
    else 
        v_theta = sqrt(mu/p)*(1+e*cos(theta));
    end

    delta_v = 2*v_theta*sin(alpha/2);

elseif delta_OM<0 && delta_i<0
    alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(abs(delta_OM)));
    cosu_i = (cos(i_f) - cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
    cosu_f = (-cos(i_i) + cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
    sinu_i = (sin(abs(delta_OM))*sin(i_f)) / sin(alpha);
    sinu_f = (sin(abs(delta_OM))*sin(i_i)) / sin(alpha);
    u_i = atan2(sinu_i, cosu_i);
    u_f = atan2(sinu_f, cosu_f);

    theta = u_i - om_i;
    om_f = u_f - theta;

     if cos(theta)>0 
        v_theta = sqrt(mu/p)*(1+e*cos(theta + pi));
    else 
        v_theta = sqrt(mu/p)*(1+e*cos(theta));
    end

    delta_v = 2*v_theta*sin(alpha/2);
end

end




