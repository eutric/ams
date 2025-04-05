function [delta_v1, delta_v2, delta_t] = bitangentTransfer(orbit_i, orbit_f, type)

a_i = orbit_i.a;
e_i = orbit_i.e;
mu = orbit_i.mu;
a_f = orbit_f.a;
e_f = orbit_f.e;

switch type 
    case 'pa'
        % pericentro -> apocentro
        rp_i = a_i*(1-e_i);
        ra_f = a_f*(1+e_f);
        rp_t = rp_i;
        ra_t = ra_f;
        a_t = (ra_t + rp_t)/2;

        vp_i = sqrt(mu)*sqrt(2/rp_i - 1/a_i);
        vp_t = sqrt(mu)*sqrt(2/rp_t - 1/a_t);
        va_f = sqrt(mu)*sqrt(2/ra_f - 1/a_f);
        va_t = sqrt(mu)*sqrt(2/ra_t - 1/a_t);

        % costo manovra 
        delta_v1 = abs(vp_t - vp_i);
        delta_v2 = abs(va_f - va_t);
        delta_t = pi*sqrt(a_t^3 / mu); % tempo di manovra 

    case 'ap'
        % apocentro -> pericentro 
        ra_i = a_i*(1+e_i);
        rp_f = a_f*(1-e_f);
        ra_t = ra_i;
        rp_t = rp_f;
        a_t = (ra_t + rp_t)/2;

        va_t = sqrt(mu)*sqrt(2/ra_t - 1/a_t);
        vp_t = sqrt(mu)*sqrt(2/rp_t - 1/a_t);
        va_i = sqrt(mu)*sqrt(2/ra_i - 1/a_i);
        vp_f = sqrt(mu)*sqrt(2/rp_f - 1/a_f);

        % costo manovra 
        delta_v1 = va_t - va_i;
        delta_v2 = vp_f - vp_t;
        delta_t = pi*sqrt(a_t^3 / mu); % tempo di manovra 

    case 'pp'
        % pericentro -> pericentro
        rp_i = a_i*(1-e_i);
        rp_f = a_f*(1-e_f);
        rp_t = rp_i;
        ra_t = rp_f;
        a_t = (ra_t + rp_t)/2;

        vp_i = sqrt(mu)*sqrt(2/rp_i - 1/a_i);
        vp_t = sqrt(mu)*sqrt(2/rp_t - 1/a_t);
        vp_f = sqrt(mu)*sqrt(2/rp_f - 1/a_f);
        va_t = sqrt(mu)*sqrt(2/ra_t - 1/a_t);

        % costo manovra 
        delta_v1 = vp_t - vp_i;
        delta_v2 = vp_f - va_t;
        delta_t = pi*sqrt(a_t^3 / mu); % tempo di manovra 

    case 'aa'
        % apocentro -> apocentro
        ra_i = a_i*(1+e_i);
        ra_f = a_f*(1+e_f);
        rp_t = ra_i;
        ra_t = ra_f;
        a_t = (ra_t + rp_t)/2;

        vp_t = sqrt(mu)*sqrt(2/rp_t - 1/a_t);
        va_t = sqrt(mu)*sqrt(2/ra_t - 1/a_t);
        va_i = sqrt(mu)*sqrt(2/ra_i - 1/a_i);
        va_f = sqrt(mu)*sqrt(2/ra_f - 1/a_f);

        % costo manovra 
        delta_v1 = vp_t - va_i;
        delta_v2 = va_f - va_t;
        delta_t = pi*sqrt(a_t^3 / mu); % tempo di monovra

end

end

