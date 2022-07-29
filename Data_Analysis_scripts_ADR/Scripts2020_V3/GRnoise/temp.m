
close all

XDATA = linspace(10,100,100);
bs = -20
xpar = [10^(bs) 0.5];
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt')
options.OptimalityTolerance =  1*10^(bs-100);
options.StepTolerance =  1*10^(bs-100);
options.FunctionTolerance =  1*10^(bs-100);
%options.Algorithm = 'trust-region-reflective';
options
%options = optimoptions(options, 'OptimalityTolerance', 10^(bs-15)); 
Model_TLS = @(x,fdata)x(1)*power(fdata,x(2)); %Model we use to fit
YDATA = Model_TLS(xpar,XDATA) + (3*10^(bs-1)).*randn(1,length(XDATA));
Y = Model_TLS(xpar,XDATA);

x0 = [(10^-(bs+2)) 1];
% Non-linear fit is done in next line!
%options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
[x,resnorm,~,exitflag,output] = lsqcurvefit(Model_TLS,x0,XDATA,YDATA,[],[],options);
f1 = figure;
ax1 = axes('XScale','linear','YScale','linear');
hold(ax1,'on')
plot(XDATA,Y,XDATA,YDATA,XDATA,Model_TLS(x,XDATA));
legend('Perfect','Generated Data','Non-Linear fit')
hold(ax1,'off')