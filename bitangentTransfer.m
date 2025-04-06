function [delta_v1, delta_v2, delta_t, orbit_bt] = bitangentTransfer(orbit_i, orbit_f, type)
% La funzione prende in ingresso parametri dell'orbita iniziale e dell'orbita finale, 
% una stringa e restituisce i due costi della manovra di trasferimento
% bitangente, il tempo di manovra e i paramentri dell'orbita di trasferimento.

% Input
% orbit_i e orbit_f sono strutture contenenti parametri caratterizzanti dell'orbita
% iniziale e finale:
%  orbit_i.a   valore del semiasse maggiore iniziale 
%  orbit_i.e   modulo del vettore eccentricità iniziale
%  orbit_f.a   valore del semiasse maggiore finale 
%  orbit_f.e   modulo del vettore eccentricità finale
%  type        stringa che indica il tipo di manovra 


% Output
% delta_v1     primo impulso della manovra
% delta_v2     secondo impulso della manovra
% delta_t      tempo di manovra
% orbit_bt     struttura contenente i parametri dell'orbita di trasferimento


% Parametri iniziali e finali 
a_i = orbit_i.a;
e_i = orbit_i.e;
mu = orbit_i.mu;
a_f = orbit_f.a;
e_f = orbit_f.e;

switch type 
    case 'pa'
        % pericentro -> apocentro
        % definizione dei raggi
        rp_i = a_i*(1-e_i);
        ra_f = a_f*(1+e_f);
        rp_t = rp_i;
        ra_t = ra_f;
        a_t = (ra_t + rp_t)/2;
        e_t = (ra_t - rp_t)/(ra_t + rp_t);

        % definizione delle velocità
        vp_i = sqrt(mu)*sqrt(2/rp_i - 1/a_i);
        vp_t = sqrt(mu)*sqrt(2/rp_t - 1/a_t);
        va_f = sqrt(mu)*sqrt(2/ra_f - 1/a_f);
        va_t = sqrt(mu)*sqrt(2/ra_t - 1/a_t);

        % costi manovra 
        delta_v1 = abs(vp_t - vp_i);
        delta_v2 = abs(va_f - va_t);
        delta_t = pi*sqrt(a_t^3 / mu); % tempo di manovra 

        % definisco orbita di trasferimento 
        orbit_bt = orbit_i;
        orbit_bt.a = a_t;
        orbit_bt.e = e_t;


    case 'ap'
        % apocentro -> pericentro 
        % definizione dei raggi
        ra_i = a_i*(1+e_i);
        rp_f = a_f*(1-e_f);
        ra_t = ra_i;
        rp_t = rp_f;
        a_t = (ra_t + rp_t)/2;
        e_t = (ra_t - rp_t)/(ra_t + rp_t);


        % definizione delle velocità
        va_t = sqrt(mu)*sqrt(2/ra_t - 1/a_t);
        vp_t = sqrt(mu)*sqrt(2/rp_t - 1/a_t);
        va_i = sqrt(mu)*sqrt(2/ra_i - 1/a_i);
        vp_f = sqrt(mu)*sqrt(2/rp_f - 1/a_f);

        % costi manovra 
        delta_v1 = va_t - va_i;
        delta_v2 = vp_f - vp_t;
        delta_t = pi*sqrt(a_t^3 / mu); % tempo di manovra 

        % definisco orbita di trasferimento 
        orbit_bt = orbit_i;
        orbit_bt.a = a_t;
        orbit_bt.e = e_t;

    case 'pp'
        % pericentro -> pericentro
        % definizione dei raggi
        rp_i = a_i*(1-e_i);
        rp_f = a_f*(1-e_f);
        rp_t = rp_i;
        ra_t = rp_f;
        a_t = (ra_t + rp_t)/2;
        e_t = (ra_t - rp_t)/(ra_t + rp_t);


        % definizione delle velocità
        vp_i = sqrt(mu)*sqrt(2/rp_i - 1/a_i);
        vp_t = sqrt(mu)*sqrt(2/rp_t - 1/a_t);
        vp_f = sqrt(mu)*sqrt(2/rp_f - 1/a_f);
        va_t = sqrt(mu)*sqrt(2/ra_t - 1/a_t);

        % costi manovra 
        delta_v1 = vp_t - vp_i;
        delta_v2 = vp_f - va_t;
        delta_t = pi*sqrt(a_t^3 / mu); % tempo di manovra 

        % definisco orbita di trasferimento 
        orbit_bt = orbit_i;
        orbit_bt.a = a_t;
        orbit_bt.e = e_t;

    case 'aa'
        % apocentro -> apocentro
        % definizione dei raggi
        ra_i = a_i*(1+e_i);
        ra_f = a_f*(1+e_f);
        rp_t = ra_i;
        ra_t = ra_f;
        a_t = (ra_t + rp_t)/2;
        e_t = (ra_t - rp_t)/(ra_t + rp_t);


        % definizione delle velocità
        vp_t = sqrt(mu)*sqrt(2/rp_t - 1/a_t);
        va_t = sqrt(mu)*sqrt(2/ra_t - 1/a_t);
        va_i = sqrt(mu)*sqrt(2/ra_i - 1/a_i);
        va_f = sqrt(mu)*sqrt(2/ra_f - 1/a_f);

        % costi manovra 
        delta_v1 = vp_t - va_i;
        delta_v2 = va_f - va_t;
        delta_t = pi*sqrt(a_t^3 / mu); % tempo di monovra

        % definisco orbita di trasferimento 
        orbit_bt = orbit_i;
        orbit_bt.a = a_t;
        orbit_bt.e = e_t;

end

end

