function delta_t = TOF(orbit, th_1, th_2)
% Time of Flight
% Algoritmo valido per orbite 0<e<1
% Problema indiretto 

a = orbit.a;
e = orbit.e;
mu = orbit.mu;
T=2*pi*sqrt(a^3/mu);




% Calcolo le anomalie eccentriche 
E_1 = 2*atan(sqrt((1-e)/(1+e)) * tan(th_1/2));
E_2 = 2*atan(sqrt((1-e)/(1+e)) * tan(th_2/2));
if th_1 > pi
    E_1 = 2*pi+E_1;
end
if th_2 > pi
    E_2 = 2*pi+E_2;
end
M_1=E_1-e*sin(E_1);
M_2=E_2-e*sin(E_2);

n= sqrt(mu/(a^3));


t_1=M_1/n
t_2=M_2/n
if th_2 < th_1
    delta_t = t_2-t_1+T;
else
    delta_t = t_2-t_1;
end
% Il tempo di volo


end
