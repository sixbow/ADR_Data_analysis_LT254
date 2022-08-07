%% Part 1: GR-noise fit. 
% This script does the main data analysis to how the interplay is between
% the TLS-noise and the GR-noise when we vary the temperature
% Author: Sietse de Boer
clc 
clear all;
close all;

%% User variables
kidn_iter = 1; %[Index] Select KID.
% IndexPopt' Here you can select which powers you want 
    
Tbath_iter = 1:14; %[Index] Number of Temperatures T_bath. Likely loop between 1:14;

% Variable (Change at your own risk!)


%% Importing Data
%
ChipInfo_path = ['..' filesep '..' filesep ]; %root path where data is, one higher than the scripts
FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D'];% This is where the
load([ChipInfo_path FFTsubsubdir filesep 'NOISE_2D.mat'])
%% Fit TLS ~2 - 20Hz

%% Fit GR-noise + TLS ~400Hz - 100KHz

%% Save these parameters in a struct/class?
% contents: C_TLS, gamma, Res_TLS_gamma
% C_GR, tau_qp, Residual_GR_Tau
% P_read, Q_i, T_bath, Q_c , KID#
% Measure for flatnes (Derivative!!)
% Or you could say that a 2 db Diffence between GR-noise and measured data.
tot_combi = length(concat_v(kidn_iter,IndexPsort,1));
f = gobjects(tot_combi,1);
ax = gobjects(tot_combi,1);

%% Plotting + save figure!
 %This has side effect that it creates a figure if it is not yet present so needs to be after you create figure.
for kidn = kidn_iter % iterate over CKIDs.
    % IndexPsort is a vector of all the indices of powers for a specific
    % KID in Increasing order. 
    %power_iter = IndexPopt(kidn); %[Index] - Popt -Number of P_read values. For now I look at 1 power. 
    power_iter = IndexPsort{kidn,1}; %[Index] - All P_read Number of P_read values. For now I look at 1 power. 
    
    for p = power_iter % This is the power iterator over a specific kid 
    f(p) = figure;
    ax(p) = axes('XScale','log','YScale','linear');
    hold(ax(p),'on')
    disp(['p = ' string(p)])
    Tcolors = colormapJetJB(length(Tbath_iter));
        for nT = Tbath_iter 
            %Fig. 1: Plotting the fractional frequency plot 
            %toplot = NOISE(p).FFTnoise{1}(:,4) > 0;
            %plot(NOISE(p).FFTnoise{1}(toplot,1),10*log10(NOISE(Pr).FFTnoise{1}(toplot,4)),...
            %         '-','color','k','LineWidth',2)
            toplot = NOISE(p).FFTnoise{nT}(:,4) > 0;
            plot(NOISE(p).FFTnoise{nT}(toplot,1),10*log10(NOISE(Pr).FFTnoise{nT}(toplot,4)),...
                    '-','color',Tcolors(nT,:),'LineWidth',1)
            xlabel('F [Hz]');ylabel('S_F/F^2 [dBc/Hz]')
            xlim([0.5,1e5]);grid on;ylim([-220,-140])
            title('')

        end % End Bath temperature.
    hold(ax(p),'off')
    export_path_graph = append('../../../Export_Figures_noGit/LT254_DA_figures/Part_1_GR/f1KID',string(kidn),'T',sprintf('%1.3f',NOISE(p).Temperature(nT)),'.png');
    exportgraphics(ax(p),export_path_graph)  
    end % End power 
end % End #KID








%% Part 2: TLS-noise coof.
% in this part of the script we want to do a TLS-noise study to find out
% how the TLS-noise depends on power. 