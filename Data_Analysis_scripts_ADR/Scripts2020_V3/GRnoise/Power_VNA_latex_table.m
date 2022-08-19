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



F0 = [4.9083 5.076 5.542 5.822 5.913 6.0052]' ;
Qc = [2.03e4 2.20e4 3.25e4 2.23e4 5.42e4 3.11e4]';
Qi = [1.43e5 8.526e4 4.94e5 5.80e5 6.1e5 6.59e5]';
Ql = [1.78e4 1.751e4 3.05e4 2.15e4 4.97e4 2.98e4]';
% The origin of this data is:S21analysis_GRV5 run on the S21 folder with th
% e power subfolder at 120mK. But getting correct indexing is not worth the
% efford right now to make this table.
% Make table of KID_N KID F0 , F0(Design) , Qc, Qi, Ql
Pindex = 5;
for kidn=1:6
    KID_N(kidn,1) = kidn;
    
    f3dbres(kidn,1) = F0(kidn,1)*(10^9)/(2*Ql(kidn,1));
end
M = [KID_N F0 Qc Qi Ql f3dbres]

csvwrite('Power_VNA_latex_table.csv',M)

function p = findp(kidn,Pindex,IndexPsort)
    p = IndexPsort{kidn,1}(Pindex);
end