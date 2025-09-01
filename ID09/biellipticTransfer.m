function [delta_v1, delta_v2, delta_v3, delta_t1, delta_t2, orbit_biel1, orbit_biel2] = biellipticTransfer(orbit_i,orbit_f, ra_t)
% Funzione che consente, tramite biellittica bitangente, di passare, in
% termini di forma, da orbit_i a orbit_f. 
% INPUT
% Name     Type    Size
% orbit_i  Struct  1x1
% orbit_f  Struct  1x1
% ra_t     Scalar  1x1
a_i=orbit_i.a;
e_i=orbit_i.e;
mu=orbit_i.mu;
a_f=orbit_f.a;
e_f=orbit_f.e;

rp_i = a_i*(1-e_i);
rp_t1= rp_i;
rp_f = a_f*(1-e_f);
rp_t2 = rp_f;
a_t1 = (ra_t + rp_t1)/2;
a_t2 = (ra_t + rp_t2)/2;

% calcolo delle velocità 
vp_i = sqrt(mu)*sqrt(2/rp_i - 1/a_i);
vp_f = sqrt(mu)*sqrt(2/rp_f - 1/a_f);
vp_t1 = sqrt(mu)*sqrt(2/rp_t1 - 1/a_t1);
vp_t2 = sqrt(mu)*sqrt(2/rp_t2 - 1/a_t2);
va_t1 = sqrt(mu)*sqrt(2/ra_t - 1/a_t1);
va_t2 = sqrt(mu)*sqrt(2/ra_t - 1/a_t2);

% calcolo il costo 
delta_v1 = abs(vp_t1 - vp_i);
delta_v2 = abs(va_t2 - va_t1);
delta_v3 = abs(vp_f - vp_t2);

% tempo di manovra 
delta_t1 = pi*sqrt(a_t1^3 / mu);
delta_t2 = pi*sqrt(a_t2^3 / mu);


% definisco i parametri delle orbite di trasferimento 
orbit_biel1 = orbit_i;
orbit_biel2 = orbit_i;
orbit_biel1.a = a_t1;
orbit_biel2.a = a_t2;
orbit_biel1.e = ra_t/a_t1 - 1;
orbit_biel2.e = ra_t/a_t2 - 1;

end

