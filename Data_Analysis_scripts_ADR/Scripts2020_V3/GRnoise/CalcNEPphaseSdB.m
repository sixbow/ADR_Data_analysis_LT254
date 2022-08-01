function NEPdark = CalcNEPphaseSdB(dFdNqp,Q,F0,Stheta,eta_pb,tau_qp,Tc_al)
%Calculates the dark NEP phase per 1 kid at a time. 1 temperature..
%multiple frequencies
% Input
% dF_{t}/dN_qp (Number)
% Q: Quality factor (Number)
% F0: res freq (Number) [Hz]
% Stheta: S_{\theta}(1,n - vector)
% tau_qp(T) (Number) [s]
% Tc : Critical temp Al (Number)[k]
% Based on: [Janssen,2014 Title: equivalence of optical... ]
% <strong>Warning</strong> the calculation of Nqp(T) is also dependent on tc in the S21
% script

dthetadNqp = (-4*Q/(F0)).*dFdNqp;
dthetadPdark = (((eta_pb)*(tau_qp))/(deltaBCS(Tc_al)))*dthetadNqp;
NEPdark = sqrt(Stheta)./dthetadPdark;
end

