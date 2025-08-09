function[deltav,deltav_r,alpha]=calcolo_velocita_man_combo(orbit_i,orbit_f,theta1,theta2)
a_i = orbit_i.a;
e_i = orbit_i.e;
i_i = orbit_i.i;
OM_i = orbit_i.OM;
om_i = orbit_i.om;
mu = orbit_i.mu;

a_f = orbit_f.a;
e_f = orbit_f.e;
i_f=orbit_f.i;
OM_f=orbit_f.OM;
om_f=orbit_f.om;

deltaOM=abs(OM_f-OM_i);

alpha=acos(cos(i_i)*cos(i_f)+sin(i_i)*sin(i_f)*cos(deltaOM));
p_i=a_i*(1-e_i^2);
p_f=a_f*(1-e_f^2);
v_r1=sqrt(mu/p_i)*e_i*sin(theta1);
v_r2=sqrt(mu/p_f)*e_f*sin(theta2);
v_theta1=sqrt(mu/p_i)*(1+e_i*cos(theta1));
v_theta2=sqrt(mu/p_f)*(1+e_f*cos(theta2));
deltav_r=abs(v_r2-v_r1);
deltav_theta=sqrt(v_theta2^2+v_theta1^2-2*v_theta1*v_theta2*cos(alpha));
deltav=sqrt(deltav_r^2+deltav_theta^2);

end

