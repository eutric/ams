function r = r_parametrica(O, theta)
%R_PARAMETRICA Summary of this function goes here
%   Detailed explanation goes here
r = O.a*(1-O.e^2)/(1+O.e*cos(theta));
end

