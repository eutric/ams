function v_r = v_radiale(O, theta)
%V_RADIALE Summary of this function goes here
%   Detailed explanation goes here
p = O.a * (1 - O.e^2);
v_r = sqrt(O.mu/p)*O.e*sin(theta);
end

