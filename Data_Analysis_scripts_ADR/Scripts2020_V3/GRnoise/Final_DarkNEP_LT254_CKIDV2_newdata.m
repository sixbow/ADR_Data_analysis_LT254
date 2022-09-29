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
%Gauconst = 1.62E-6;%8E-7; % const before Q*alpha_k/V    V is in um^3 % Beuninggg..
Gauconst = 8.2E-7;%8E-7; %
AlVol = [66 66 66 132 132 132];% Detector volume in um^3 
%tau_qp() = We do this in the code
kidn_iter = 1:6;
kidn_iter2_2_2 = 1:3;
kidn_iter4_4_4 = 4:6;

nT_iter = 1:14;
%</Input>
freq_vec = real(CrossPSDNOISE(1).CrossPSD{1,1}(:,1));
NEP_Matrix = zeros(length(kidn_iter),length(nT_iter),length(freq_vec));
NEP_Matrix_max = zeros(length(kidn_iter),length(nT_iter),length(freq_vec));
NEP_Matrix_min = zeros(length(kidn_iter),length(nT_iter),length(freq_vec));
dthetadPdark = zeros(length(kidn_iter),length(nT_iter),length(freq_vec));
dthetadNqp = zeros(length(kidn_iter),length(nT_iter),length(freq_vec));
NEP_GR_Theory_Matrix = zeros(length(kidn_iter),length(nT_iter));
for kidn = kidn_iter2_2_2 % Loop over KIDS
    for nT = nT_iter % Loop over temperatures
    Stheta_current_lin = dbtolin(real(CrossPSDNOISE(kidn).CrossPSD{1,nT}(:,3)));%Stheta is always real..
    %Stheta_current_lin = dbtolin(real(CrossPSDNOISE(kidn).CrossPSD{1,nT}(:,3))./10);%Stheta is always real..
    
    tau_qp_current = CrossPSDFit(kidn).tau(nT);
    tau_qp_current_max = CrossPSDFit(kidn).taumax(nT);
    tau_qp_current_min = CrossPSDFit(kidn).taumin(nT);
    [NEP_Matrix(kidn,nT,:),dthetadPdark(kidn,nT,:),dthetadNqp(kidn,nT,:)] = CalcNEPphaseSdB(dthetaindNqp2_2_2(kidn),Q(kidn),F0(kidn),Stheta_current_lin,freq_vec,eta_pb,tau_qp_current,Tc_al); % 3 dimention is freq.
    [NEP_Matrix_max(kidn,nT,:),~,~] = CalcNEPphaseSdB(dthetaindNqp2_2_2(kidn),Q(kidn),F0(kidn),Stheta_current_lin,freq_vec,eta_pb,tau_qp_current_max,Tc_al);
    [NEP_Matrix_min(kidn,nT,:),~,~] = CalcNEPphaseSdB(dthetaindNqp2_2_2(kidn),Q(kidn),F0(kidn),Stheta_current_lin,freq_vec,eta_pb,tau_qp_current_min,Tc_al);
    
    alpha_k(kidn) = dthetaindNqp2_2_2(kidn).*AlVol(kidn)./(Gauconst.*Q(kidn));
    alpha_k_high(kidn) = 1.62*alpha_k(kidn);
    alpha_k_low(kidn) = 0.617*alpha_k(kidn);
    %GR limited NEP
    % Possible to add NEP gr
    [~,~,~,NEP_GR_Theory_Matrix(kidn,nT)] = getNqp_Tbath(Tc_al,eta_pb,AlVol(kidn),NOISE(kidn).Temperature(nT));
    end
end

for kidn = kidn_iter4_4_4 % Loop over KIDS
    for nT = nT_iter % Loop over temperatures
    Stheta_current_lin = dbtolin(real(CrossPSDNOISE(kidn).CrossPSD{1,nT}(:,3)));%Stheta is always real..
    %Stheta_current_lin = dbtolin(real(CrossPSDNOISE(kidn).CrossPSD{1,nT}(:,3))./10);%Stheta is always real..
    
    tau_qp_current = CrossPSDFit(kidn).tau(nT);
    [NEP_Matrix(kidn,nT,:),dthetadPdark(kidn,nT,:),dthetadNqp(kidn,nT,:)] = CalcNEPphaseSdB(dthetaindNqp4_4_4(kidn),Q(kidn),F0(kidn),Stheta_current_lin,freq_vec,eta_pb,tau_qp_current,Tc_al); % 3 dimention is freq.
    [NEP_Matrix_max(kidn,nT,:),~,~] = CalcNEPphaseSdB(dthetaindNqp4_4_4(kidn),Q(kidn),F0(kidn),Stheta_current_lin,freq_vec,eta_pb,tau_qp_current_max,Tc_al);
    [NEP_Matrix_min(kidn,nT,:),~,~] = CalcNEPphaseSdB(dthetaindNqp4_4_4(kidn),Q(kidn),F0(kidn),Stheta_current_lin,freq_vec,eta_pb,tau_qp_current_min,Tc_al);
    
    alpha_k(kidn) = dthetaindNqp4_4_4(kidn).*AlVol(kidn)./(Gauconst*Q(kidn));
    alpha_k_high(kidn) = 1.62*alpha_k(kidn);
    alpha_k_low(kidn) = 0.617*alpha_k(kidn);
    [~,~,~,NEP_GR_Theory_Matrix(kidn,nT)] = getNqp_Tbath(Tc_al,eta_pb,AlVol(kidn),NOISE(kidn).Temperature(nT));
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
        %plot(freq_vec,squeeze(NEP_Matrix_max(kidn,nT,:)),':','Color',plt_color_map(nT,:),'LineWidth',1);
        %plot(freq_vec,squeeze(NEP_Matrix_min(kidn,nT,:)),':','Color',plt_color_map(nT,:),'LineWidth',1);
        if nT>8
        yline(NEP_GR_Theory_Matrix(kidn,nT),'--','Color',plt_color_map(nT,:),'LineWidth',2,'HandleVisibility','off')
        end
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
kidn = 2;
ax2 = axes('XScale','log','YScale','log');
hold(ax2, 'on');

for nT = nT_iter
    plot(freq_vec,squeeze(NEP_Matrix(kidn,nT,:)),'Color',plt_color_map(nT,:),'LineWidth',2);
    %plot(freq_vec,squeeze(NEP_Matrix_max(kidn,nT,:)),':','Color',plt_color_map(nT,:),'LineWidth',1);
    %plot(freq_vec,squeeze(NEP_Matrix_min(kidn,nT,:)),':','Color',plt_color_map(nT,:),'LineWidth',1);
    if nT>8
    yline(NEP_GR_Theory_Matrix(kidn,nT),'--','Color',plt_color_map(nT,:),'LineWidth',2,'HandleVisibility','off')
    end
end
hold(ax2, 'off');
title(append('CKID #',string(kidn)));
ylabel('$NEP_{dark}$ ($W/\sqrt{Hz}$)','Interpreter','latex');
xlabel('f [Hz]','Interpreter','latex');
for nT=nT_iter 
    legendTvalues{nT} = sprintf('%1.3f',NOISE(kidn).Temperature(nT));
end
hTl2(kidn) = legend(legendTvalues,'location','eastOutside');
title('Dark NEP as function of $T_{bath} | G3C8 $','Interpreter','latex')
title(hTl2(kidn),'T_{bath} [K]')
hTl2(kidn).ItemTokenSize = [3,10];
grid on

f3 = figure
hold on
bar(1:6,alpha_k)
er = errorbar(1:6,alpha_k,alpha_k_low-alpha_k,alpha_k_high-alpha_k,'HandleVisibility','off')
er.Color = [0 0 0]
plot(1:6,[0.32 0.32 0.32 0.21 0.21 0.21],'o','LineWidth',1.5)
alphak_noise = [0.317 0.317 0.317 0.1758 0.1758 0.1758];
alphak_noise_high = [0.356 0.356 0.356 0.1973 0.1973 0.1973];
alphak_noise_low = [0.283 0.283 0.283 0.1567 0.1567 0.1567];
plot(1:6,alphak_noise,'x','LineWidth',1.5,'Color','r')
er2 = errorbar(1:6,alphak_noise,alphak_noise_low-alphak_noise,alphak_noise_high-alphak_noise,'HandleVisibility','off')
er2.Color = [1 0 0]
legend('\alpha_{k,dark}','\alpha_{k,design}','\alpha_{k,noise}')
er.LineStyle = 'none';
title('\alpha_{k} from dark responsivity')
ylabel('\alpha_{k}')
xlabel('CKID #')
hold off




%exportgraphics(f1,'../../../Export_Figures_noGit/Final_DarkNEP_images/f1NEP_dark.png')
%exportgraphics(f1,'../../../Export_Figures_noGit/Final_DarkNEP_images/f1NEP_dark.pdf')
%exportgraphics(f2,'../../../Export_Figures_noGit/Final_DarkNEP_images/f2NEP_dark.png')
%exportgraphics(f2,'../../../Export_Figures_noGit/Final_DarkNEP_images/f2NEP_dark.pdf')
%exportgraphics(f3,'../../../Export_Figures_noGit/Final_DarkNEP_images/f3NEP_dark.png')
%exportgraphics(f3,'../../../Export_Figures_noGit/Final_DarkNEP_images/f3NEP_dark.pdf')