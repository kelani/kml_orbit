clc
clear all
close all 
si = 2;
mu = 3;
x=-10:0.2:20;
gauss = (1/(si*sqrt(2*pi)))*(exp(-((x-mu).^2)/(2*(si^2))));
plot (gauss);
hold on;
si = 30;
mu = 40;
x=-10:0.2:20;
plot (gauss);
hold off;
