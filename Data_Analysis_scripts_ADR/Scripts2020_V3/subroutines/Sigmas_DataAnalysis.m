function [sigma1,sigma2,ds1dnqp,ds2dnqp, dxdn_alpha_1] = Sigmas_DataAnalysis(frequ, N0, Delta0, T)
% uses analytical expressions to calcualte sigma1/sigman, sigma 2/sigman (2.16, 2.17 Pieter Thesis) 
% all expressions only valid for Delta = Dealta0!
% copied from Jochem's scriptset
% IN
% frequ angular frequency of resonator
% N0 is ensity of sates at fermi level in eV-1um-3 (as in the main script
% Deltain is in (J), used for sigma.
% T temperature value T (K)
% OUT
% sigma 1 and sigma2 are both normalized to sigman (as in 2.16, 2.17)
%  dsigma1dnqp, dsigma2dnqp are both normalized to sigman and per um^3
% (2.18 and 2.19 Pieter thesis)
%  nqp is #quasiparticles/um^3!!!

%NB: 
e_c=1.602e-19;          %single electron charge
kb=1.3806503e-23;       %boltzmann contstant J/K
hbar=6.60607e-34/2/pi;  %1.054571e-34;  %reduced plancks constant

Beta = 2;

NOJ = N0/e_c;  %now in /J/um^3

%(2.16, 2.17 Pieter Thesis) 
ksi=hbar.*frequ./(2*kb.*T);
sigma1=4*Delta0./(hbar.*frequ).*exp(-Delta0./(kb.*T)).*sinh(ksi).*besselk(0,ksi);
sigma2=pi*Delta0./(hbar.*frequ).*( 1-2.*exp(-Delta0./(kb.*T)).*exp(-ksi).*besseli(0,ksi)  );

%(2.18 and 2.19 Pieter thesis)    
ds1dnqp=1/NOJ./(hbar*frequ).*sqrt(2*Delta0/pi/kb./T).*sinh(ksi).*besselk(0,ksi);
ds2dnqp=-1*pi/2/NOJ./(hbar*frequ).*(1+sqrt(2*Delta0/pi/kb./T).*exp(-ksi).*besseli(0,ksi));

%
abssigma=1*abs(sigma1-1i*sigma2);%normalized to sigman
dxdn_alpha_1    = -1*Beta*0.25/abssigma * ds2dnqp; %Pieter 2.27+3.38, frequ. response per qp for alpha = 1
end

