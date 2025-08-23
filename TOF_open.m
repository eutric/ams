function delta_t = TOF_open(orbit, th_1, th_2)
% Time of Flight
% Algoritmo valido per orbite 0<e<1
% Problema indiretto 

a = orbit.a;
e = orbit.e;
mu = orbit.mu;




% Calcolo le anomalie eccentriche 
F_1 = 2*atanh(sqrt((e-1)/(1+e)) * tan(th_1/2));
F_2 = 2*atanh(sqrt((e-1)/(1+e)) * tan(th_2/2));
if th_1 > pi
    F_1 = 2*pi+F_1; % E_1 sarà negativo, a th=pi varrebbe E=-pi, aumentando
    % E_1 diminuisce ed assume il valore negativo a partire da 0 (quindi va
    % bene)
end
if th_2 > pi
    F_2 = 2*pi+F_2;
end

% M_1=F_1-e*sin(F_1);
% M_2=F_2-e*sin(F_2);

n=sqrt(-mu/(a^3));

t1=1/n*(e*sinh(F_1)-F_1);
t2=1/n*(e*sinh(F_2)-F_2);

delta_t=t2-t1;

% t_1=M_1/n;
% t_2=M_2/n;
% if th_2 < th_1
%     delta_t = t_2-t_1+T;
% else
%     delta_t = t_2-t_1;
% end


end


