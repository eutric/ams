
clear all
clc
close all
teta=linspace(0,10,10);
g=plot(0,0,'ro');
for i=1:length(teta)
    set(g,"YData",teta(i))
    pause (1)
end