clc
clear all
close all

XDATA = logspace(1,5,100);
bs = -7;
xpar = [10^(bs) -0.5];
%options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt')
%options = optimoptions('lsqcurvefit', 'TolFun', 1e-80, 'TolX', 1e-80, 'Display',   'iter', 'DiffMinChange', .01);
%options.MaxFunctionEvaluations = 10000;
%options.MaxIterations = 10000;
%options.OptimalityTolerance =  1*10^(bs-100);
%options.StepTolerance =  1*10^(bs-100);
%options.FunctionTolerance =  1*10^(bs-100);
%options.Algorithm = 'trust-region-reflective';
%options = optimoptions(options, 'OptimalityTolerance', 10^(bs-15)); 
Model_TLS_lin = @(x,fdata)x(1)*power(fdata,x(2)); %Model we use to fit
YDATA = (Model_TLS_lin(xpar,XDATA) + (3*10^(bs-2)).*randn(1,length(XDATA)));
Y = Model_TLS_lin(xpar,XDATA);

x0 = [10^(-7) 1];

x = nonlinfitSdB(XDATA,YDATA,Model_TLS_lin,x0);

f1 = figure;
ax1 = axes('XScale','linear','YScale','linear');
hold(ax1,'on')
Y_reconstructed = Model_TLS_lin(x,XDATA);
plot(XDATA,Y,XDATA,YDATA);
plot(XDATA,Y_reconstructed);
legend('Perfect','Generated Data','Non-Linear fit');
hold(ax1,'off')