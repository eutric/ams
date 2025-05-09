clear all
clc
close all
%definito k come rapporto ra1/ra2 
mu=398600;
O_start.a=24400.0;
O_start.e=0;
O_start.i=0.104700;
O_start.OM=1.514000;
O_start.om=3.107000;
O_start.mu=mu;
th_start=1.665000;
[rr,vv]=par2car(O_start,th_start);
h=norm(cross(rr,vv));
%caratterizzo le orbite con h ed e scelti a partire dai dati prime
P=h^2/mu;

O_1=O_start;
rp1=P/(1+O_start.e);
ra1=P/(1-O_start.e);
O_1.a=(ra1+rp1)/2;
O_1.e=sqrt(1-P/O_1.a);
plotOrbit(O_1,0,2*pi,0.01,'-')
hold on
N=100;
cost=zeros(N,2);
ii=0;
k=linspace(1,40,N);
for kk=k
    ii=ii+1;
O_2=O_1;
ra2=ra1*kk;
% O_2.e=-P/ra2+1;
% rp2=P/(1+O_2.e); 
rp2=kk*rp1;
O_2.a=(ra2+rp2)/2;
plotOrbit(O_2,0,2*pi,0.01,'--')
[delta_v1, delta_v2, delta_t, orbit_bt,th0,thf] = bitangentTransfer(O_1, O_2, 'pa');
cost(ii,1)=delta_v1+delta_v2;
cost_be=0;
cost_be_vect=zeros(100,N);
min=0;
jj=0;
for j=linspace(ra2,100*ra2,100)
    jj=jj+1;
[delta_v1, delta_v2, delta_v3, delta_t1, delta_t2] = biellipticTransfer(O_1,O_2, j);
cost_be=delta_v1+delta_v2+delta_v3;
if cost_be<min || min==0
    min=cost_be;
end
cost_be_vect(jj,ii)=cost_be;


end
cost(ii,2)=min;


end
figure
plot(k,cost(:,1));
hold on
plot(k,cost(:,2));
legend('costo bitangente','costo migliore biellittica')
figure

plot(cost_be_vect)
hold on
