%generazione orbita
clear all
close all
clc

e=0.3;
h=3;
mu=3;
P=orbit(e,mu,h);
plot(P(1,:),P(2,:))
hold on
plot(0,0,'o',Color='k',LineWidth=1.5)
grid on
%% legge oraria su orbita
clear all
close all
clc
%genero orbita
e=0.3;
h=3;
mu=3;
p=h^2/mu;
a=p/(1-e^2);
rp=p/(1+e);
P=orbit(e,mu,h);
k=plot(P(1,:),P(2,:));
hold on
b=plot(0,0,'o',Color='k',LineWidth=1.5);
grid on
% tempo
T=linspace(0,10,1000);
N=length(T);
n=0.5; %rad/s
tetVect=zeros(N,1);
tet0=0;
E0=0;
period=2*pi*sqrt(a^3/mu);
count=1;
for i=2:N
    count=count+1;
    t=T(i);
    if abs(tetVect(i-1))<pi
        fE=@(Ev)n*t-Ev+e*sin(Ev);
        E=fzero(fE,E0);
        E0=E;
        ftet=@(te)tan(E/2)-sqrt((1-e)/(1+e))*tan(te/2);
        tet=fzero(ftet,tet0);
        tetVect(i)=tet;
        tet0=tet;
    else 
        break
    end
end
teta_2=2*pi*ones(length(tetVect(1:count-2)),1)-flip(tetVect(1:count-2));
teta=[tetVect(1:count-2);teta_2];
teta=[teta;teta];
%provo plot
Z=plot(rp,0,'ro');

for i=2:length(teta)
    r=p/(1+e*cos(teta(i)));
    rx=r*cos(teta(i));
    ry=r*sin(teta(i));
    set(Z,"XData",rx,"YData",ry);
    pause(0.01)
% end
    % rxvect=[rxvect;rx];
    % ryvect=[ryvect;ry];
end
% for i=1:length(ryvect)
%     






