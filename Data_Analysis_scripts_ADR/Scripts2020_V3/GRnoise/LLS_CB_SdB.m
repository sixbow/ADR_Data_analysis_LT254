function coof = LLS_CB_SdB(XDATA,YDATA,func,coof0,n_num)
%nonlinfitSdB <strong>nonlinfitSdB(XDATA,YDATA,func,coof0)</strong>
% This version normalizes the YDATA to make it stable because it uses trust
% region standard. This needs to be adapted for the problem at hand. 

YDATA_N = YDATA./n_num;
coof0(1) = (coof0(1)./n_num);
[x,~,~,~,~] = lsqcurvefit(func,coof0,XDATA,YDATA_N);%[],[],options);
coof(1) = x(1).*n_num ; %This makes use of that th
coof(2)= x(2);
end

