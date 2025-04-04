function delta_t = TOF(orbit, th_1, th_2)
% Time of Flight
% Algoritmo valido per orbite 0<e<1
% Problema indiretto 

a = orbit.a;
e = orbit.e;
mu = orbit.mu;
T=2*pi*sqrt(a^3/mu);

if th_1 <= pi && th_2>pi
    





% Calcolo le anomalie eccentriche 
E_1 = 2*atan(sqrt((1-e)/(1+e)) * tan(th_1/2));
E_2 = 2*atan(sqrt((1-e)/(1+e)) * tan(th_2/2));

M_1=E_1-e*sin(E_1);
M_2=E_2-e*sin(E_2);

n= sqrt(a^3 /mu);


t_1=M_1/n;
t_2=M_2/n;

% Il tempo di volo

if th_1 > th_2
    
    delta_t=t_2-t_1+T;

   
