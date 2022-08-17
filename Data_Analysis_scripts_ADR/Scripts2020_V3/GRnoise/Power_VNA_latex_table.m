% Simple script to generate the table of resonance frequencies extracted
% out of NOISE

clc
clear all 
close all 

ChipInfo_path = ['..' filesep '..' filesep ]; %root path where data is, one higher than the scripts       
FFTsubdir = ['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_120mK' filesep 'FFT' filesep 'Power']; % 120mK    %
filename = 'Noise_P.mat';
NOISE = load([ChipInfo_path FFTsubdir filesep filename],'NOISE').NOISE;
IndexPsort = load([ChipInfo_path FFTsubdir filesep filename],'IndexPsort').IndexPsort;


% Make table of KID_N KID F0 , F0(Design) , Qc, Qi, Ql
Pindex = 5;
for kidn=1:6
    KID_N(kidn,1) = kidn;
    F0(kidn,1) = NOISE(findp(kidn,Pindex,IndexPsort)).Fres;
    Qc(kidn,1) = NOISE(findp(kidn,Pindex,IndexPsort)).Qc;
    Qi(kidn,1) = NOISE(findp(kidn,Pindex,IndexPsort)).Qi;
    Ql(kidn,1) = NOISE(findp(kidn,Pindex,IndexPsort)).Ql;
end
M = [KID_N F0 Qc Qi Ql]

csvwrite('Power_VNA_latex_table.csv',M)

function p = findp(kidn,Pindex,IndexPsort)
    p = IndexPsort{kidn,1}(Pindex);
end