function v_t = v_theta(O, theta)
%V_THETA Summary of this function goes here
%   Detailed explanation goes here
p = O.a*(1-O.e^2);
v_t = sqrt(O.mu/p)*(1+O.e*cos(theta));
end

