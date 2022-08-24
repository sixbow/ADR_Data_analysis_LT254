function [nqp, tauqp, Nqp, NEP] = getNqp_Tbath(Tc,eta_pb,V,Tbath,tau0)
%calculates the nqp, tau_qp for given Tbath
%
% would need an update to give Cref and Cro as input to deakl with
% different alu types
% INPUT
% Tc of alu
% eta_pb is the efficiency of qp creation. 0.4 for alu
% V in um^3
% 
% OUTPUT
% nqp in qp/um^3, 
% tau = lifetime in sec, 
% N qp is tota # qp
% NEP = NEP for nqp, tau qp and V. If V not given, it takes V=104um^3, for
% LT020 Chip 3
%
%
% constants
e_c=1.602e-19;          %single electron charge
N0 = 1.72e10/e_c;             %density of states in /J/um^3
kb=1.3806503e-23;       %boltzmann contstant J/K
%%%%%%%%%%%% fixed parameters alumninium %%%%%%%%%%%%
Delta_J=1.76*kb*Tc;  %in J. 
%%%%%%%%%%%%%%%

nqp = 2*N0 * (2*pi*kb*Tbath*Delta_J).^0.5 .* exp(-Delta_J./(kb*Tbath));

if nargin == 4
    tau0=458e-9;
end

Cref = 0.5*tau0 * kb * Tc * N0 * 1/(1.76)^2;

tauqp   = Cref./nqp;

Nqp     = nqp*V;
NEP     = 2*Delta_J/eta_pb * (nqp*V/tauqp).^0.5;

end