% This code is used as a subroutine that changes
% CrossPSDNoise.mat to remove the offset of the TLS noise in the signal to
% improve fitting
% TODO: Rename  CrossPSDNoise_2D -> CrossPSDNoise_2D_INPUT
% Step 1: importing the mat file 
% Step 2: getting the offset of the TLS-noise
% Step 3: Subtract and save
% (Original file is not touched!) 
% Author: Sietse de Boer 
clc
clear all;%This should be commented out for performance reasons.
close all;
ChipInfo_path = ['..' filesep '..' ]; %root path where data is, one higher than the scripts
FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D_Popt'];% This is where the
Outputfolder = 'CPSDMinusTLS';
Outputfolderdir = [ChipInfo_path,filesep,FFTsubsubdir,filesep,Outputfolder] ; 
% CrossPSDNoise_2D file is. %FFTsubdir = [filesep 'Noise_Powers_165mK' filesep 'FFT' filesep 'Power'];     %
% Loading files!
matfile = 'Noise_2D.mat';
matfile2 = 'CrossPSDNoise_2D';
matfile3 = 'CrossPSDFit_2D';
load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','IndexP_sub_opt','KIDnumbers');
load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile2],'CrossPSDNOISE');
CPSDfolder = [ChipInfo_path, filesep, 'CrossPSD2D'] ;
if ~exist(Outputfolderdir, 'dir')
       mkdir(Outputfolderdir)
end % This makes a subfolder called CPSDMinusTLS which will be were the files for this script will be created
% will be making for loop, But for now 
%++++INPUTS++++++++++
% options = optimoptions('lsqcurvefit', 'TolFun', 1e-80, 'TolX', 1e-80, 'Display',   'iter', 'DiffMinChange', .01);
% options.MaxFunctionEvaluations = 1000;
% options.MaxIterations = 1000;
%++++Switches to choose what to compensate. 
SW1_subtract_oneoverf = 1;
SW2_ring_mix_noise = 1; % Choose 1 for subtraction and 0 -> do nothing. Choose -1 for addition(Why would you want that ?? U mad?)
SW = [SW1_subtract_oneoverf SW2_ring_mix_noise];
%---/Switches

dfg = 125000; % [Hz]
fring = 5.5E9/(pi*18000); %Is really dependent on T, But lets put it like this for now!
begin_data_point = 3;
end_data_point = 30; % We only want to fit the first part of the TLS noise since there it is dominant
%figshow =[{'on'} {'on'} {'on'}] %determines which figures to show!
figshow =[{'off'} {'off'} {'off'} {'on'}] %determines which figures to show!
%figshow =[{'off'} {'off'} {'off'}] %determines which figures to show!
%---\INPUTS----------
%% Analysis
close all
for kidn = 1 %Number of KIDS(1:6)
    for nT =1:2:14%Number of Temperatures. Likely loop between 1:14  ~length():
        p = 1; % Popt only has 1 power so for now we do just one. 
        %Data
        C_freq = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1);
        C_S_CPSD_Re = real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2));
        C_S_CPSD_Im = imag(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2));
        yb = mean(C_S_CPSD_Re(121:131));
        %Fig. 1: Shows the complete CPSD with all components. Saves them
        %also!
        f1 = figure;
        ax1 = axes('XScale','log','YScale','log');
        hold(ax1,'on')
        plot(C_freq,C_S_CPSD_Re,'-o','MarkerFaceColor','r','Color','r');
        plot(C_freq,-C_S_CPSD_Re,'-o','MarkerFaceColor','g','Color','g');
        plot(C_freq,C_S_CPSD_Im,'-o','MarkerFaceColor', 'b','Color','b');
        plot(C_freq,-C_S_CPSD_Im,'-o','MarkerFaceColor', 'y','Color','y');
        legend('Re{CPSD}','-Re{CPSD}','Im{CPSD}','-Im{CPSD}')
        hold(ax1,'off')
        grid on
        set(f1, 'Visible', figshow{1});
        title(append('KID',string(kidn),'T',sprintf('%1.3f',NOISE(p).Temperature(nT))))
        export_path_graph = append('../../../Export_Figures_noGit/TLS_surpression/f0FiTSubKID',string(kidn),'T',sprintf('%1.3f',NOISE(p).Temperature(nT)),'.png');
        exportgraphics(ax1,export_path_graph)
        
        %Ring-time detuning mixing noise compensation.
        %Fig. 2: Shows the effect of compensation of the ring time on the
        %real axis.
        f2 = figure;
        ax2 = axes('XScale','log','YScale','log');
        hold(ax2,'on')
        plot(C_freq,C_S_CPSD_Re,'-o','MarkerFaceColor','r','Color','r');
        plot(C_freq,-C_S_CPSD_Re,'-o','MarkerFaceColor','g','Color','g');
        plot(C_freq,ring_mix_noise(yb,dfg,C_freq,fring),'--');
        S_model_ring = ring_mix_noise(yb,dfg,C_freq,fring);
        C_S_CPSD_Re_corr_ring= C_S_CPSD_Re-S_model_ring;
        plot(C_freq,C_S_CPSD_Re_corr_ring,'-x','Color','r','LineWidth',4);
        plot(C_freq,-C_S_CPSD_Re_corr_ring,'-x','Color','g','LineWidth',4);
        hold(ax2,'off')
        grid on
        set(f2,'Visible',figshow{2});
        
        %Fitting the slope using a nonlinear fit S(f)= C_{TLS}\frac{1}{f^{a}}
        %Status: Showing that it works!
        mean_interval = mean(C_S_CPSD_Re(begin_data_point:end_data_point))
        x0 = [mean_interval];
        Model_oneoverf = @(x,fdata)x(1)*power(fdata,-1); %Model we use to fit
        xcoof_Re = nonlinfitSdB(C_freq(begin_data_point:end_data_point),C_S_CPSD_Re(begin_data_point:end_data_point),Model_oneoverf,x0);
        S_model_oneoverf = Model_oneoverf(xcoof_Re,C_freq);
        C_S_CPSD_Re_corr_oneoverf = C_S_CPSD_Re-S_model_oneoverf;
        
        %Fig. 3: Shows how the power law fits in the data, Shows the
        %effect of compensation.
        f3 = figure;
        ax3 = axes('XScale','log','YScale','log');
        hold(ax3,'on')
        plot(C_freq,C_S_CPSD_Re,'-o','MarkerFaceColor','r','Color','r');
        plot(C_freq,-C_S_CPSD_Re,'-o','MarkerFaceColor','g','Color','g');
        plot(C_freq,Model_oneoverf(xcoof_Re,C_freq),'Color','r','LineStyle','--');
        plot(C_freq,-Model_oneoverf(xcoof_Re,C_freq),'Color','r','LineStyle','--');
        plot(C_freq,C_S_CPSD_Re_corr_oneoverf,'Color','magenta','LineStyle','--','Marker','x');
        plot(C_freq,-C_S_CPSD_Re_corr_oneoverf,'Color','cyan','LineStyle','--','Marker','x');
        xline(C_freq(begin_data_point),'--')
        xline(C_freq(end_data_point),'--')
        grid on
        hold(ax3,'off')
        grid on
        set(f3, 'Visible',figshow{3});
        title(append('KID',string(kidn),'T',sprintf('%1.3f',NOISE(p).Temperature(nT))))
        export_path_graph = append('../../../Export_Figures_noGit/TLS_surpression/f3FiTSubKID',string(kidn),'T',sprintf('%1.3f',NOISE(p).Temperature(nT)),'.png');
        exportgraphics(ax3,export_path_graph)
           
        %<----RESULT----> Resulting compensated function
        C_PSD_Re_corr = C_S_CPSD_Re-SW(1).*S_model_oneoverf -SW(2).*S_model_ring;
        CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2) = C_PSD_Re_corr + 1i.*C_S_CPSD_Im;
        %</----RESULT----> Resulting compensated function
        
        %Fig. 3: Shows how the power law fits in the data, Shows the
        %effect of compensation.
        f4 = figure;
        ax4 = axes('XScale','log','YScale','log');
        hold(ax4,'on')
        plot(C_freq,C_S_CPSD_Re,'-o','MarkerFaceColor','r','Color','r');
        plot(C_freq,-C_S_CPSD_Re,'-o','MarkerFaceColor','g','Color','g');
        plot(C_freq,SW(1).*S_model_oneoverf +SW(2).*S_model_ring,'Color','r','LineStyle','--');
        plot(C_freq,-SW(1).*S_model_oneoverf -SW(2).*S_model_ring,'Color','g','LineStyle','--');
        %plot(C_freq,S_model_oneoverf,'Color','r','LineStyle','--');
        %plot(C_freq,-S_model_oneoverf,'Color','g','LineStyle','--');
        plot(C_freq,C_PSD_Re_corr,'Color','magenta','LineStyle','--','Marker','x');
        plot(C_freq,-C_PSD_Re_corr,'Color','cyan','LineStyle','--','Marker','x');
        xline(C_freq(begin_data_point),'--')
        xline(C_freq(end_data_point),'--')
        grid on
        hold(ax4,'off')
        grid on
        set(f4, 'Visible',figshow{4});
        title(append('KID',string(kidn),'T',sprintf('%1.3f',NOISE(p).Temperature(nT))))
        export_path_graph = append('../../../Export_Figures_noGit/TLS_surpression/f4FiTSubKID',string(kidn),'T',sprintf('%1.3f',NOISE(p).Temperature(nT)),'.png');
        exportgraphics(ax4,export_path_graph)   
    end
end

%% Closing
save([Outputfolderdir,filesep,matfile],'NOISE','IndexP_sub_opt','KIDnumbers');
save([Outputfolderdir,filesep,matfile2],'CrossPSDNOISE');
disp('FitandSubtract: Jobs done!');
% disp('Press any key to quit and clean up!')
% pause
% clear all