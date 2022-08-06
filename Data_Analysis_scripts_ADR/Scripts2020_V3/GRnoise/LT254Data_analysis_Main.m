%% Part 1: GR-noise fit. 
% This script does the main data analysis to how the interplay is between
% the TLS-noise and the GR-noise when we vary the temperature
% Author: Sietse de Boer
clc 
clear vars;
close all;

%% Importing Data
%
ChipInfo_path = ['..' filesep '..' ]; %root path where data is, one higher than the scripts
FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D'];% This is where the
load([ChipInfo_path FFTsubsubdir filesep 'NOISE_2D.mat'])
%% Fit TLS ~2 - 20Hz

%% Fit GR-noise + TLS ~400Hz - 100KHz

%% Save these parameters in a struct/class?
% contents: C_TLS, gamma, Res_TLS_gamma
% C_GR, tau_qp, Residual_GR_Tau
% P_read, Q_i, T_bath, Q_c , KID#
% Measure for flatnes (Derivative!!)


%% Plotting + save figure!



%% Part 2: TLS-noise coof.
% in this part of the script we want to do a TLS-noise study to find out
% how the TLS-noise depends on power. 