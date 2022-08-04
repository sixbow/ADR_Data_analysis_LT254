function [S,debug] = Create3D_S_Gau(CrossPSDNOISE,IndexP_sub_opt,kidn,nT,p)
%Building the S matrix of the [Gau,Mazin paper 2007 APL 90,102507] 
% INPUT: CROSSPSDNOISE(Struct),kidn(number - index MKID),nT(number - index Temp),p(number - index Readout power)
% OUT: S(2x2xfreq) = [S_II S_IQ ; conj(S_IQ) S_QQ]

% Building the matrices in the third dimention.
SII = permute(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,3),[3 2 1]);
%disp(SII);
debug = SII;
SIQ = permute(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2),[3 2 1]);
SQQ = permute(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,4),[3 2 1]);
% Constructing the matrix!
S = [SII,SIQ;conj(SIQ),SQQ];   
end

