function S = ring_mix_noise(yb,dfg,f,fring)
%Function described in appendix of Pieter de Visser dissertation. About
%mixing noise due to detuning of the generator freq.
%S = yb.*(1+((dfg.*(dfg+f))./fring))./(2.*(1+((dfg+f)./fring).^2)); Not
%working
S = yb.*0.25.*abs(gamma_PdV(dfg,f,fring)-conj(gamma_PdV(dfg,-f,fring)));
end

