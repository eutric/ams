function delta_t = TOF(orbit, th_1, th_2)
% Time of Flight
% Algoritmo valido per orbite 0<e<1
% Problema indiretto 

a = orbit.a;
e = orbit.e;
mu = orbit.mu;


% Calcolo le anomalie eccentriche 
E_1 = 2*atan(sqrt((1-e)/(1+e)) * tan(th_1/2));
E_2 = 2*atan(sqrt((1-e)/(1+e)) * tan(th_2/2));


% Il tempo di volo
if th_2 > th_1
    delta_t = sqrt(a^3 /mu)* ((E_2 - E_1) - e*(sin(E_2) - sin(E_1)));
else 
    delta_t = sqrt(a^3 /mu)* ((E_2 - E_1) - e*(sin(E_2) - sin(E_1))) + 2*pi*sqrt(a^3/mu);

end

end
