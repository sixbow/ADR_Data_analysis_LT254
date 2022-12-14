%% runs all scripts can take time..
clc
clear all 
close all
%% Part 1. To find P_opt Comment in line 17 , Comment out line 18
% After this has been found put Popt files in cd .. -> 2D_Popt folder.
% So now the analysis proceeds on only Popt files.
run('NOISEanalysis_PdepV7.m')
clear all 
close all
%% (Done) Does PSD analysis for all temperatures based on the time domain data that is present
% Also calculates the crossPSD.
% Known issue. Could not open file ? Makegoodfigure...
run('NOISEanalysis_2DV1.m')
%%
%pause; % If you want to save some figures uncomment this line.
clear all 
close all
%% (Done) Does CPSD analysis for all temperatures. (Runtime: ~10min)
% Fixed:Known issue. Could not open file ? Makegoodfigure...
% This makes the CPSD
run('TDanalysisV2.m')
% Where are figures saved?
%pause; % If you want to save some figures uncomment this line.
%%
clear all 
close all
%% (Done) Remove TLS component - Attempt at removing the TLS by fitting and subtracting.
run('FitandSubtractTLS_CPSD_SDBV1.m')
% This saves in seperate outputfolder:/CPSDMinusTLS
%% (Doing) Fitting using Jochem script
% use the outputfolder of previous script as input.:/CPSDMinusTLS
% Outputs Figures are in: ..\..\CrossPSD2D\KID6P_87_crossPSD_Tdep
run('CrossPSDJBV2.m')
%%
clear all 
close all
%% Part 2 (Done)S21 analysis (F0(T),Qi(T))
% Added S21input.txt to timedomain folder , tab seperated KID#  F0[GHz] Q
% L_alu[um]
run('S21analysis_GRV5.m')



%% Main analysis script
clear all 
close all 
run('Final_DarkNEP_LT254_CKID.m')


disp('Finished CKID analysis!')