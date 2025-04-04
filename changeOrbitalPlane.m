function [delta_v,om_f, teta,O_end] = changeOrbitalPlane(orbit_i, orbit_f)
% funzione che cambia il piano di manovra 
a=orbit_i.a;
e=orbit_i.e;
i_i = orbit_i.i;
OM_i = orbit_i.OM;
om_i = orbit_i.om;
mu=orbit_i.mu;
i_f=orbit_f.i;
OM_f=orbit_f.OM;

delta_OM = OM_f - OM_i;
delta_i = i_f - i_i;
p = a*(1-e^2);
%creo orbita d arrivo
O_end=orbit_i;
O_end.OM=OM_f;
O_end.i=i_f;
if delta_OM>0 && delta_i>0
    alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(abs(delta_OM)));
    cosu_i = (-cos(i_f) + cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
    cosu_f = (cos(i_i) - cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
    sinu_i = (sin(abs(delta_OM))*sin(i_f)) / sin(alpha);
    sinu_f = (sin(abs(delta_OM))*sin(i_i)) / sin(alpha);
    u_i = atan2(sinu_i, cosu_i);
    u_f = atan2(sinu_f, cosu_f);

    teta = u_i - om_i;
    om_f = u_f - teta;

    if cos(teta)>0 
    v_teta = sqrt(mu/p)*(1+e*cos(teta + pi));
    else 
        v_teta = sqrt(mu/p)*(1+e*cos(teta));
    end

    delta_v = 2*v_teta*sin(alpha/2);


elseif delta_OM>0 && delta_i<0
    alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(abs(delta_OM)));
    cosu_i = (cos(i_f) - cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
    cosu_f = (-cos(i_i) + cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
    sinu_i = (sin(abs(delta_OM))*sin(i_f)) / sin(alpha);
    sinu_f = (sin(abs(delta_OM))*sin(i_i)) / sin(alpha);
    u_i = atan2(sinu_i, cosu_i);
    u_f = atan2(sinu_f, cosu_f);

    teta = 2*pi - u_i - om_i;
    om_f = 2*pi - u_f - teta;

    if cos(teta)>0 
    v_teta = sqrt(mu/p)*(1+e*cos(teta + pi));
    else 
        v_teta = sqrt(mu/p)*(1+e*cos(teta));
    end

    delta_v = 2*v_teta*sin(alpha/2);

elseif delta_OM<0 && delta_i>0
    alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(abs(delta_OM)));
    cosu_i = (-cos(i_f) + cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
    cosu_f = (cos(i_i) - cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
    sinu_i = (sin(abs(delta_OM))*sin(i_f)) / sin(alpha);
    sinu_f = (sin(abs(delta_OM))*sin(i_i)) / sin(alpha);
    u_i = atan2(sinu_i, cosu_i);
    u_f = atan2(sinu_f, cosu_f);

    teta = 2*pi - u_i - om_i;
    om_f = 2*pi - u_f - teta;
     if cos(teta)>0 
    v_teta = sqrt(mu/p)*(1+e*cos(teta + pi));
    else 
        v_teta = sqrt(mu/p)*(1+e*cos(teta));
    end

    delta_v = 2*v_teta*sin(alpha/2);

elseif delta_OM<0 && delta_i<0
    alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(abs(delta_OM)));
    cosu_i = (cos(i_f) - cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
    cosu_f = (-cos(i_i) + cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
    sinu_i = (sin(abs(delta_OM))*sin(i_f)) / sin(alpha);
    sinu_f = (sin(abs(delta_OM))*sin(i_i)) / sin(alpha);
    u_i = atan2(sinu_i, cosu_i);
    u_f = atan2(sinu_f, cosu_f);

    teta = u_i - om_i;
    om_f = u_f - teta;

     if cos(teta)>0 
    v_teta = sqrt(mu/p)*(1+e*cos(teta + pi));
    else 
        v_teta = sqrt(mu/p)*(1+e*cos(teta));
    end

    delta_v = 2*v_teta*sin(alpha/2);
end

end




