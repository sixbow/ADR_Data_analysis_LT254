% Over here I tie the analysis together and plot the NEP of the CKID
% devices.
clc
close all 
clearvars
load('../../Data_LT254_Sietse/LT254_Sietse_Chip11/Noise_vs_T/FFT/2D_Popt/CPSDMinusTLS/CrossPSDFit_2D.mat')
load('../../Data_LT254_Sietse/LT254_Sietse_Chip11/Noise_vs_T/FFT/2D_Popt/CrossPSDNoise_2D.mat')
% Load spectra...
S21_path = '../../Data_LT254_Sietse/LT254_Sietse_Chip11/S21/2D_Popt/all_Tdep.csv';
S21_results = table2array(readtable(S21_path , "NumHeaderLines",1));
dFdNqp =  S21_results(:,11);%Warning! Why is this negative valued!
%<Input>
Q = S21_results(:,3); % vec for all KIDS
F0 = S21_results(:,3); % for all kids
eta_pb = 0.4;
Tc_al = 1.182;
%tau_qp() = We do this in the code
kidn_iter = 1:6;
nT_iter = 1:14;
%</Input>
freq_vec = real(CrossPSDNOISE(1).CrossPSD{1,1}(:,1));
NEP_Matrix = zeros(length(kidn_iter),length(nT_iter),length(freq_vec));
for kidn = kidn_iter % Loop over KIDS
    for nT = nT_iter % Loop over temperatures
    Stheta_current_lin = dbtolin(real(CrossPSDNOISE(kidn).CrossPSD{1,nT}(:,3))./10);%Stheta is always real..
    tau_qp_current = CrossPSDFit(kidn).tau(nT);
    NEP_Matrix(kidn,nT,:) = CalcNEPphaseSdB(dFdNqp(kidn),Q(kidn),F0(kidn),Stheta_current_lin,eta_pb,tau_qp_current,Tc_al); % 3 dimention is freq.
    end
end
% Fig 1. 
f1 = figure;
for kidn=kidn_iter
    ax(kidn) = subplot(2,3,kidn,'XScale','log','YScale','log')
    hold(ax(kidn), 'on');
    %color = colormapJetJB(length(IndexP_sub_opt{kidn}));
    for nT = nT_iter
        plot(freq_vec,-squeeze(NEP_Matrix(kidn,nT,:)));
    end
    hold(ax(kidn), 'off');
    title(['CKID #' string(kidn)])
    grid on
end


