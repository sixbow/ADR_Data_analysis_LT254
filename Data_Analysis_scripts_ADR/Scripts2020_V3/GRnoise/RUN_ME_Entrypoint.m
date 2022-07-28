%% runs all scripts can take time..
clc
clear all 
close all
%% To find P_opt Comment in line 17 , Comment out line 18
% After this has been found put Popt files in cd .. -> 2D_Popt folder.
% So now the analysis proceeds on only Popt files.
run('NOISEanalysis_PdepV7.m')
clear all 
close all
%%(Done) Does PSD analysis for all temperatures. 
% Known issue. Could not open file ? Makegoodfigure...
run('NOISEanalysis_2DV1.m')
%pause; % If you want to save some figures uncomment this line.
clear all 
close all
%%(Done) Does PSD analysis for all temperatures. 
% Known issue. Could not open file ? Makegoodfigure...
run('TDanalysisV2.m')
%pause; % If you want to save some figures uncomment this line.
clear all 
close all
%% Remove TLS component - Attempt at removing the TLS by fitting and subtracting.
run('FitandSubtractTLS_CPSD_SDBV1.m')

%% Fitting using Jochem script
run('CrossPSDJBV2.m')
clear all 
close all

%% S21 analysis (F0(T),Qi(T))


%% Main analysis script



disp('Finished CKID analysis!')