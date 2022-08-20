% Over here I tie the analysis together and plot the NEP of the CKID
% devices.
clc
close all 
clearvars
load('../../Data_LT254_Sietse/LT254_Sietse_Chip11/Noise_vs_T/FFT/2D_Popt/CPSDMinusTLS/CrossPSDFit_2D.mat')
load('../../Data_LT254_Sietse/LT254_Sietse_Chip11/Noise_vs_T/FFT/2D_Popt/CrossPSDNoise_2D.mat')
load('../../Data_LT254_Sietse/LT254_Sietse_Chip11/Noise_vs_T/FFT/2D_Popt/NOISE_2D.mat','NOISE');
% Load spectra...
S21_path = '../../Data_LT254_Sietse/LT254_Sietse_Chip11/S21/2D_Popt/all_Tdep2_2_2.csv';
S21_results2_2_2 = table2array(readtable(S21_path , "NumHeaderLines",1));
S21_path = '../../Data_LT254_Sietse/LT254_Sietse_Chip11/S21/2D_Popt/all_Tdep4_4_4.csv';
S21_results4_4_4 = table2array(readtable(S21_path , "NumHeaderLines",1));
dthetaindNqp2_2_2 =  S21_results2_2_2(:,11);%Warning! Why is this negative valued!(Sietse: i think this is normal!)
dthetaindNqp4_4_4 =  S21_results4_4_4(:,11);%Warning! Why is this negative valued!(Sietse: i think this is normal!)
%<Input>
Q = S21_results2_2_2(:,3); % vec for all KIDS
F0 = S21_results2_2_2(:,6); % for all kids
eta_pb = 0.4;
Tc_al = 1.182;
Gauconst = 8E-7; % const before Q*alpha_k/V    V is in um^3 % Beuninggg..
AlVol = [66 66 66 132 132 132];% Detector volume in um^3 
%tau_qp() = We do this in the code
kidn_iter = 1:6;
kidn_iter2_2_2 = 1:3;
kidn_iter4_4_4 = 4:6;

nT_iter = 1:14;
%</Input>
freq_vec = real(CrossPSDNOISE(1).CrossPSD{1,1}(:,1));
NEP_Matrix = zeros(length(kidn_iter),length(nT_iter),length(freq_vec));
dthetadPdark = zeros(length(kidn_iter),length(nT_iter),length(freq_vec));
dthetadNqp = zeros(length(kidn_iter),length(nT_iter),length(freq_vec));
for kidn = kidn_iter2_2_2 % Loop over KIDS
    for nT = nT_iter % Loop over temperatures
    Stheta_current_lin = dbtolin(real(CrossPSDNOISE(kidn).CrossPSD{1,nT}(:,3)));%Stheta is always real..
    %Stheta_current_lin = dbtolin(real(CrossPSDNOISE(kidn).CrossPSD{1,nT}(:,3))./10);%Stheta is always real..
    
    tau_qp_current = CrossPSDFit(kidn).tau(nT);
    [NEP_Matrix(kidn,nT,:),dthetadPdark(kidn,nT,:),dthetadNqp(kidn,nT,:)] = CalcNEPphaseSdB(dthetaindNqp2_2_2(kidn),Q(kidn),F0(kidn),Stheta_current_lin,freq_vec,eta_pb,tau_qp_current,Tc_al); % 3 dimention is freq.
    alpha_k(kidn) = dthetaindNqp2_2_2(kidn).*AlVol(kidn)./(Gauconst*Q(kidn));
    end
end

for kidn = kidn_iter4_4_4 % Loop over KIDS
    for nT = nT_iter % Loop over temperatures
    Stheta_current_lin = dbtolin(real(CrossPSDNOISE(kidn).CrossPSD{1,nT}(:,3)));%Stheta is always real..
    %Stheta_current_lin = dbtolin(real(CrossPSDNOISE(kidn).CrossPSD{1,nT}(:,3))./10);%Stheta is always real..
    
    tau_qp_current = CrossPSDFit(kidn).tau(nT);
    [NEP_Matrix(kidn,nT,:),dthetadPdark(kidn,nT,:),dthetadNqp(kidn,nT,:)] = CalcNEPphaseSdB(dthetaindNqp4_4_4(kidn),Q(kidn),F0(kidn),Stheta_current_lin,freq_vec,eta_pb,tau_qp_current,Tc_al); % 3 dimention is freq.
    alpha_k(kidn) = dthetaindNqp4_4_4(kidn).*AlVol(kidn)./(Gauconst*Q(kidn));
    end
end
% Fig 1. 
f1 = figure('units','normalized','outerposition',[0 0 1 1]);
plt_color_map = colormapJetJB(length(nT_iter)); %Source Jochem code.
for kidn=kidn_iter
    
   
    ax(kidn) = subplot(2,3,kidn,'XScale','log','YScale','log');
    hold(ax(kidn), 'on');
    
    for nT = nT_iter
        plot(freq_vec,squeeze(NEP_Matrix(kidn,nT,:)),'Color',plt_color_map(nT,:),'LineWidth',2);
    end
    hold(ax(kidn), 'off');
    title(append('CKID #',string(kidn)));
    ylabel('$NEP_{dark}$ ($W/\sqrt{Hz}$)','Interpreter','latex');
    xlabel('f [Hz]','Interpreter','latex');
    for nT=nT_iter 
        legendTvalues{nT} = sprintf('%1.3f',NOISE(kidn).Temperature(nT));
    end
    hTl(kidn) = legend(legendTvalues,'location','eastOutside');
    title(hTl(kidn),'T_{bath} [K]')
    hTl(kidn).ItemTokenSize = [3,10];
    grid on
end
f2 = figure('units','normalized','outerposition',[0 0 1 1]);
kidn = 1;
ax2 = axes('XScale','log','YScale','log');
hold(ax2, 'on');

for nT = nT_iter
    plot(freq_vec,squeeze(NEP_Matrix(kidn,nT,:)),'Color',plt_color_map(nT,:),'LineWidth',2);
end
hold(ax2, 'off');
title(append('CKID #',string(kidn)));
ylabel('$NEP_{dark}$ ($W/\sqrt{Hz}$)','Interpreter','latex');
xlabel('f [Hz]','Interpreter','latex');
for nT=nT_iter 
    legendTvalues{nT} = sprintf('%1.3f',NOISE(kidn).Temperature(nT));
end
hTl2(kidn) = legend(legendTvalues,'location','eastOutside');
title('Dark NEP as function of $T_{bath}$','Interpreter','latex')
title(hTl2(kidn),'T_{bath} [K]')
hTl2(kidn).ItemTokenSize = [3,10];
grid on

f3 = figure
bar(alpha_k)






exportgraphics(f1,'../../../Export_Figures_noGit/Final_DarkNEP_images/f1NEP_dark.png')
exportgraphics(f1,'../../../Export_Figures_noGit/Final_DarkNEP_images/f1NEP_dark.pdf')
exportgraphics(f1,'../../../Export_Figures_noGit/Final_DarkNEP_images/f2NEP_dark.png')
exportgraphics(f1,'../../../Export_Figures_noGit/Final_DarkNEP_images/f2NEP_dark.pdf')