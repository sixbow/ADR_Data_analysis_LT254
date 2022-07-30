function output = gamma_PdV(dfg,f,fring)
%<strong>Gamma</strong> function as described in the dissertation of Pieter de Visser. 
% input parameters:
% dfg = detuning freq relative to resonance frequency of MKID [Hz]
% f = modulation freq
% fring = F0/pi*Q , Resonator ring time
output = (1+1i.*(dfg./fring))./(1+(1i.*((dfg+f)./fring)));
end

