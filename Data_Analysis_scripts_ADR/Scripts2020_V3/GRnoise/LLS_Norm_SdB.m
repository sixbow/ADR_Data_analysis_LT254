function coof = LLS_CB_SdB(XDATA,YDATA,func,coof0,lb,ub,normalization_vector)
%nonlinfitSdB <strong>nonlinfitSdB(XDATA,YDATA,func,coof0)</strong>
% This version normalizes the YDATA to make it stable because it uses trust
% region standard. This needs to be adapted for the problem at hand. 

n_num = mean(YDATA);
YDATA_N = YDATA./n_num;
coof0(1) = (coof0(1)./n_num);
lb(1) = lb(1)/n_num;
ub(1) = ub(1)/n_num;% Normalizing everything..
[x,~,~,~,~] = lsqcurvefit(func,coof0,XDATA,YDATA_N,lb,ub);%[],[],options);
coof(1) = x(1).*n_num ;
coof(2)= x(2);
end
