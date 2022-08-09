klier % clearing..
N = 100;
a= 2;
b = 4;
c = 0.5;
x = linspace(-10,10,N);
f = @(x)myfunc1(x,a,b);
g = @(x)myfunc2(x,c);

fming = @(x)f(x)-g(x);
% World
hold on
plot(x,f(x),x,g(x));
grid on
% finding the intersect!
[i1(1),i1(2)] = findintersect_SdB(f,g,-2);
[i2(1),i2(2)] = findintersect_SdB(f,g,10);
plot(i1(1),i1(2),'x','LineWidth',2)
plot(i2(1),i2(2),'x','LineWidth',2)



%
%%

function out = myfunc1(x,a,b)
out = a.*x + b;
end 
function out = myfunc2(x,c)
out = c.*(x.^2);
end