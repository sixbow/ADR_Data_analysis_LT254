function coof = nonlinfitSdB(XDATA,YDATA,func,coof0)
%nonlinfitSdB <strong>nonlinfitSdB(XDATA,YDATA,func,coof0)</strong>
% This version normalizes the YDATA to make it stable.
n_num = mean(YDATA);
YDATA_N = YDATA./n_num;
coof0_N = [(coof0(1)./n_num) coof0(2)];
[x,~,~,~,~] = lsqcurvefit(func,coof0_N,XDATA,YDATA_N);%[],[],options);
coof = [x(1).*n_num x(2)];
end

