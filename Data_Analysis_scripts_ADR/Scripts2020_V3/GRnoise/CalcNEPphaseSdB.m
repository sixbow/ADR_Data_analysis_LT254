function [NEPdark,dthetadPdark,dthetadNqp] = CalcNEPphaseSdB(dFdNqp,Q,F0,Stheta,f,eta_pb,tau_qp,Tc_al)
%Calculates the dark NEP phase per 1 kid at a time. 1 temperature..
%multiple frequencies
% Input
% dF_{t}/dN_qp (Number)
% Q: Quality factor (Number)
% F0: res freq (Number) [Hz]
% f: modulation frequency [Hz]
% Stheta: S_{\theta}(1,n - vector)
% tau_qp(T) (Number) [s]
% Tc : Critical temp Al (Number)[k]
% Based on: [Janssen,2014 Title: equivalence of optical... ]
% <strong>Warning</strong> the calculation of Nqp(T) is also dependent on tc in the S21
% script

dthetadNqp = (-4*Q/(F0)).*dFdNqp;
sprintf('Q = %1.3e',Q)
sprintf('F0 = %1.3e',F0)
dthetadPdark = (((eta_pb)*(tau_qp))/(deltaBCS(Tc_al)))*dthetadNqp;
sprintf('dthetadPdark = %1.3e',dthetadPdark)
NEPdark = (sqrt(Stheta)./dthetadPdark).*sqrt(1+(2.*pi.*tau_qp.*f).^2);

end

