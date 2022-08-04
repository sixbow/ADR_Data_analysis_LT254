function [O,E,va,vb,Saa,Sbb] = Sdiag_Gau(CrossPSDNOISE,IndexP_sub_opt,kidn,nT,p)
%diagonalizing the real part of the S matrix as shown in[Gau,Mazin paper 2007 APL 90,102507] 
% INPUT: CROSSPSDNOISE(Struct),kidn(number - index MKID),nT(number - index Temp),p(number - index Readout power)
% OUT: O (Orthogonal rotation matrix -> O'=inv(O) ), E(diagonal matrix of eigenvalues)
% such that real(S)= O*E*O'

ReS = real(Create3D_S_Gau(CrossPSDNOISE,IndexP_sub_opt,kidn,nT,p));
% Here we use the Create3D_S_Gau function to construct the matrix
size_S = size(ReS); 
for i=1:size_S(3)
    [O(:,:,i),E(:,:,i)] = eig(ReS(:,:,i));
    va(:,i) = squeeze(O(:,1));
    vb(:,i) = squeeze(O(:,2));
    Saa(:,:,i) = E(1,1);
    Sbb(:,:,i) = E(2,2);
end
end

