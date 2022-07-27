function [NEP,nqp,Nqp] = NEPquick2(tauqp_exp,V,tau0_sp)
% [NEP,tauqp_exp,Nqp] = NEPquick(nqp,V)
% n = qpparticles / um^3, 10-100
% tauqp_exp = experimental lifetime in sec
% tau0_sp   = kaplan, single particle related tau0. 
% Old results from Pieter
% are for 2 particles, so tau0_sp = tau0_PdV*2;
% checked with PdV chapteer 6 data
% corrected for wrong factor 2

if nargin == 0
    tauqp_exp = 0.5e-3;
    V = 1000;
    tau0_PdV    = 458e-9;           	%from PdV p120, we accidently related this to the experimntal tau and tau_exp x 2 = tau_single particle. same for tau_0
    tau0_sp = tau0_PdV*2;
elseif nargin == 2
    tau0_PdV    = 458e-9;           	%from PdV p120, we accidently related this to the experimntal tau and tau_exp x 2 = tau_single particle. same for tau_0
    tau0_sp     = tau0_PdV*2;
   
end


%constants
kb      = 1.3806503e-23;    %boltzmann contstant J/K
e_c     = 1.602e-19;          %single electron charge

%alu properties
N0      = 1.72*10^10/e_c;       %particles um-3 J-1
Tc      = 1.25;
eta_pb  = 0.4;

tauqp_sp = 2*tauqp_exp;     %correct single particle lifetime is 2x the measured one!
%resulting params
Cref    = tau0_sp * kb * Tc * N0 * 1/(2*1.76)^2; %sec/um^3; relation between tau and nqp. replaced tauo with tau0_sp

Delta   = 1.76*kb*Tc;
    
nqp   = Cref./tauqp_sp; %this is INDEPENDENT the factor 2 story due to single particle parameters vs experiemntal 2 particle parameters as they compensate
Nqp     = nqp*V;
NEP     = 2*Delta/eta_pb * sqrt(Nqp/tauqp_exp);