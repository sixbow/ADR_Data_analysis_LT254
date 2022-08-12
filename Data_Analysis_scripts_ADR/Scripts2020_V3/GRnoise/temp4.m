clc
clear all
close all
N = 100;
XDATA = linspace(1,10,N);
a = 2;
f1 = figure;
ax = axes();
for i=1:10
YDATA = XDATA.*a + (3.*randn(1,N) -1.5 ) - 4 ;
coof0 = 1.45;

f = @(b,xdata)myfunc(xdata,a,b);

[Cout,~,~,~,~] = lsqcurvefit(f,coof0,XDATA,YDATA);


hold(ax,'on')
plot(XDATA,YDATA,XDATA,f(Cout,XDATA));
hold(ax,'off')
end
% Functions end!!
function out = myfunc(b,xdata,a,b)
out = a.*xdata + b;
end
