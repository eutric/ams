function [P] = orbit(e,mu,h)
teta=[0:1:360];
P=zeros(2,361);
for t=1:361
    r=(h.^2/mu)/(1+e.*cos(teta(t)*pi/180));
    P(1,t)=r*cos(teta(t)*pi/180);
    P(2,t)=r*sin(teta(t)*pi/180);
end

