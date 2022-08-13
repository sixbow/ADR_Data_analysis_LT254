clc
clear all;
close all;

f = @(x,C)C(1).*x+C(2);
g = @(x,D)D.*x.^2;
f1 = @(x)f(x,[2 3]);
g1 = @(x)g(x,0.3);
[x0,y0] = findintersect_SdB2(f1,g1,0.2);
hold on
plot(linspace(-10,10,100),f1(linspace(-10,10,100)))
plot(linspace(-10,10,100),g1(linspace(-10,10,100)))
plot(x0,y0,'x')
