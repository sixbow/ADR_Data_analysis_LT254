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
options = optimoptions('lsqcurvefit', 'TolFun', 1e-80, 'TolX', 1e-80, 'Display',   'iter', 'DiffMinChange', .01);
options.MaxFunctionEvaluations = 1000;
options.MaxIterations = 1000;
%figshow =[{'on'} {'on'} {'on'}] %determines which figures to show!
figshow =[{'off'} {'off'} {'off'}] %determines which figures to show!
%figshow =[{'off'} {'off'} {'off'}] %determines which figures to show!
%% Analysis
close all
for kidn = 1:6 %Number of KIDS(1:6)
    for nT =1:14%Number of Temperatures. Likely loop between 1:14  ~length():
        p = 1; % Popt only has 1 power so for now we do just one. 
        %Data
        C_freq = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1);
        C_S_CPSD_Re = real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2));
        C_S_CPSD_Im = imag(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2));
        
        f1 = figure;
        ax1 = axes('XScale','log','YScale','log');
        hold(ax1,'on')
        plot(C_freq,C_S_CPSD_Re,'-o','MarkerFaceColor','r','Color','r');
        plot(C_freq,-C_S_CPSD_Re,'-o','MarkerFaceColor','g','Color','g');
        plot(C_freq,C_S_CPSD_Im,'-o','MarkerFaceColor', 'b','Color','b');
        plot(C_freq,-C_S_CPSD_Im,'-o','MarkerFaceColor', 'y','Color','y');
        legend('Re{CPSD}','-Re{CPSD}','Im{CPSD}','-Im{CPSD}')
        begin_data_point = 1; 
        end_data_point = 100; % We only want to fit the first part of the TLS noise since there it is dominant
        hold(ax1,'off')
        grid on
        set(f1, 'Visible', figshow{1});
        title(append('KID',string(kidn),'T',sprintf('%1.3f',NOISE(p).Temperature(nT))))
        export_path_graph = append('../../../Export_Figures_noGit/TLS_surpression/f0FiTSubKID',string(kidn),'T',sprintf('%1.3f',NOISE(p).Temperature(nT)),'.png');
        exportgraphics(ax1,export_path_graph)
        
        
        f2 = figure;
        ax2 = axes('XScale','log','YScale','log');
        hold(ax2,'on')
        yb = 1E-10;
        dfg = 125000; % [Hz]
        fring = 5.5E9/(pi*20000); %Is really dependent on T, But lets put it like this for now!
        plot(C_freq,C_S_CPSD_Re,'-o','MarkerFaceColor','r','Color','r');
        plot(C_freq,-C_S_CPSD_Re,'-o','MarkerFaceColor','g','Color','g');
        plot(C_freq,ring_mix_noise(yb,dfg,C_freq,fring),'--');
        C_S_CPSD_Re_corr= C_S_CPSD_Re-ring_mix_noise(yb,dfg,C_freq,fring);
        plot(C_freq,C_S_CPSD_Re_corr,'-x','Color','r','LineWidth',4);
        plot(C_freq,-C_S_CPSD_Re_corr,'-x','Color','g','LineWidth',4);
        hold(ax2,'off')
        grid on
        set(f2, 'Visible',figshow{2});
        
        
        f3 = figure;
        ax3 = axes('XScale','log','YScale','log');
        hold(ax3,'on')
        %Insert plots
        hold(ax3,'off')
        grid on
        set(f3, 'Visible',figshow{3});
        
        
        
        
%         f1 = figure
%         % -- Not used currently
%         %x0 = [(10^-8) 1];
%         %Model_TLS = @(x,fdata)x(1)*power(fdata,-x(2)); %Model we use to fit
%         % Non-linear fit is done in next line!
%         %options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
%         %[x,resnorm,~,exitflag,output] = lsqcurvefit(Model_TLS,x0,C_freq(begin_data_point:end_data_point),C_S_CPSD_Re(begin_data_point:end_data_point),[],[],options);
%         % Plotting non-linear least squares approach..
%         %plot(C_freq,Model_TLS(x,C_freq));
%         % -- \Not used currently
%         xline(C_freq(begin_data_point),'--')
%         xline(C_freq(end_data_point),'--')
%         grid on
%         % Trying how the Linear goes with fitting in log-log! 
%         % Real part
%         Xlog_fitr = log10(C_freq(begin_data_point:end_data_point)); % f = e^x
%         Ylog_fitr_Re = log10(C_S_CPSD_Re(begin_data_point:end_data_point));
%         log_fit_Re = polyfit(Xlog_fitr,Ylog_fitr_Re,1);
%         % Imag Part
%         Xlog_fitr = log10(C_freq(begin_data_point:end_data_point)); % f = e^x
%         Ylog_fitr_Im = log10(C_S_CPSD_Im(begin_data_point:end_data_point));
%         log_fit_Im = polyfit(Xlog_fitr,Ylog_fitr_Im,1);
%         LinearModel = @(coof,f)power(10,coof(2)).*power(f,coof(1));
%         LogSpaceModel = @(coof,x) coof(1).*x + coof(2);
%         %Plotting the Linear in LogLog
%         plot(C_freq,LinearModel(log_fit_Re,C_freq));
%         plot(C_freq,LinearModel(log_fit_Im,C_freq));
%         TLS_corrected_S_CPSD_Re = C_S_CPSD_Re - LinearModel(log_fit_Re,C_freq); %+ ones(1,length(Current_freq))*min(Current_S_CPSD)  ;
%         TLS_corrected_S_CPSD_Im = C_S_CPSD_Im; %+ ones(1,length(Current_freq))*min(Current_S_CPSD)  ;
%         plot(C_freq,TLS_corrected_S_CPSD_Re);
%         plot(C_freq,TLS_corrected_S_CPSD_Im);
%         legend('Cross-PSD','Begin fit','End fit','Lin in loglogspace')
%         title(sprintf('KID# %d T= %1.3f K P_{read} %d dBm ',kidn,NOISE(p).Temperature(nT),NOISE(kidn).ReadPower))
%         hold(ax1,'off')
%         set(f1, 'Visible', 'off');
%         export_path_graph = append('../../../Export_Figures_noGit/FiTSubKID',string(kidn),'T',sprintf('%1.3f',NOISE(p).Temperature(nT)),'.png');
%         exportgraphics(ax1,export_path_graph)
%         
% %         f2a = figure;
% %         ax2a = axes('XScale','linear','YScale','linear');
% %         hold(ax2a,'on')
% %         plot(Xlog_fitr_Re,Ylog_fitr_Re);
% %         plot(Xlog_fitr_Re,LogSpaceModel(log_fit,Xlog_fitr_Re))
% %         %plot(Current_freq(begin_data_point:end_data_point),LinearModel(log_fit,Current_freq(begin_data_point:end_data_point)));
% %         hold(ax2a,'off')
% %         set(f2a, 'Visible', 'off');
% %         grid on
% 
%         f2b = figure;
%         ax2b = axes('XScale','linear','YScale','linear');
%         hold(ax2b,'on')
%         plot(C_freq(begin_data_point:end_data_point),C_S_CPSD_Re(begin_data_point:end_data_point));
%         title('Linear fitting range Real part CPSD')
%         hold(ax2b,'off')
%         set(f2b, 'Visible', 'off');
%         
%         f2 = figure;
%         ax2 = axes('XScale','linear','YScale','linear');
%         hold(ax2,'on')
%         plot(C_freq,C_S_CPSD_Im);
%         title('Linear Im part CPSD')
%         hold(ax2,'off')
%         set(f2, 'Visible', 'on');
% 
%         f3 = figure;
%         ax3 = axes('XScale','linear','YScale','linear');
%         hold(ax3,'on')
%         plot(C_freq(begin_data_point:end_data_point),Current_S_CPSD(begin_data_point:end_data_point));
%         plot(C_freq(begin_data_point:end_data_point),Model_TLS(x,C_freq(begin_data_point:end_data_point)));
%         title('What is happening in linear space in the fitting window')
%         hold(ax3,'off')
%         set(f3, 'Visible', 'off');
%         
%         %plot(Current_freq,Model_TLS(Current_freq));
%         f4 = figure;
%         ax4 = axes('XScale','linear','YScale','linear');
%         hold(ax4,'on')
%         hold(ax4,'off')
%         set(f4, 'Visible', 'off');
%         
          CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2) = C_S_CPSD_Re_corr + 1i.*C_S_CPSD_Im;
        
    end
end
%%


%% Closing

save([Outputfolderdir,filesep,matfile],'NOISE','IndexP_sub_opt','KIDnumbers');
save([Outputfolderdir,filesep,matfile2],'CrossPSDNOISE');
disp('Jobs done!');
% disp('Press any key to quit and clean up!')
% pause
% clear all