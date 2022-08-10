%% Part 1: GR-noise fit. 
% This script does the main data analysis to how the interplay is between
% the TLS-noise and the GR-noise when we vary the temperature
% Author: Sietse de Boer
klier % Clearing...
% Importing Data
%
addpath('..\subroutines')
ChipInfo_path = ['..' filesep '..' filesep ]; %root path where data is, one higher than the scripts
FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D'];% This is where the
load([ChipInfo_path FFTsubsubdir filesep 'NOISE_2D.mat'])
freq = NOISE(1).FFTnoise{1,1}(:,1);
% User variables

kidn_iter = 2; %[Index] Select KID.
% IndexPopt' Here you can select which powers you want 
    
Tbath_iter = 1:3:14; %[Index] Number of Temperatures T_bath. Likely loop between 1:14;

% Variable (Change at your own risk!)

% TLS fit range 
TLSfitmin = 0.8;% [Hz]
TLSfitmax = 30;% [Hz]
[~,TLSf_min_i] = min(abs(freq-TLSfitmin)); % This gives the min index
[~,TLSf_max_i] = min(abs(freq-TLSfitmax)); % This gives the max index
C_v_TLS_0 = [(10^-15) 1.5];% CTLS at 1 Hz , power with what  the power law decreases. TLS ~0.5
C_v0_dBlog = [-150 1.5];% CTLS at 1 Hz , power with what  the power law decreases. TLS ~0.5
%TLS_lb = [10^-25 0.1];
%TLS_ub = [10^-8 2];
% Combined fit range: GR model + TLS noise model = approx GR model at high
% freq
CBfitmin = 1000;% [Hz]
CBfitmax = 150000;%  [Hz]
[~,CBf_min_i] = min(abs(freq-CBfitmin)); % "
[~,CBf_max_i] = min(abs(freq-CBfitmax)); % "
C_v_CB_0 = [(10^-17) 0.0001];% CTLS at 1 Hz , power with what  the power law decreases. TLS ~0.5

fguess = 2000; % [Hz] initial guess for the intersection point of the TLS and GR noise. 

% CB_lb = [10^-20 0.000000000000001];
% CB_ub = [10^-15 1];
CB_lb = [10^-20 (1/(2*pi*CBfitmax))];
CB_ub = [10^-15 (1/(2*pi*CBfitmin))];
%
tot_combi = length(concat_v(kidn_iter,IndexPsort,1));
f = gobjects(tot_combi,1);
ax = gobjects(tot_combi,1);



% Define fitting model 
% Here I will define the models that we will be using to fit the data. 
% This we will use to obtain a non-linear fit of the data. This gives heavy
% weight on the parts of the graph that are high on log log plot since here
% the difference is the biggest. 
Model_oneoverf = @(C_v,fdata)C_v(1).*power(fdata,-1*C_v(2)); %Model we use to fit C_v is the constants vector that we are trying to find.
Model_GR = @(C_v,fdata) C_v(1)./((1+power((2.*pi.*fdata.*C_v(2)),2).*(1))); %Model we use to fit C_v is the constants vector that we are trying to find.
%Model_CB = @(C_v,fdata) C_TLS_corr.*power(fdata,-1*gamma_TLS)+C_v(1)./(1+power((2.*pi.*fdata.*C_v(2)),2)); %Model we use to fit C_v is the constants vector that we are trying to find.

%These models below are meant to be used in log space. Maybe later f = 10^x
% --> Later..
Model_oneoverf_dBlog = @(C_v,xdata)10*log10(C_v(1)) -10.*C_v(2).*xdata; %Model we use to fit C_v is the constants vector that we are trying to find.
%Model_GR_logspace = @(C_v,xdata) 10*log10(C(1))+10.*log10(1./(1+power((2.*pi.*power(10,x).*C(2)),2))); %Model we use to fit C_v is the constants vector that we are trying to find.

%Linear helper line
Model_line = @(C_v,xdata)C_v(1).*xdata+C_v(2);

%% Part 1: Sec 1 | Plotting + save figure!
% REMEMBER to load 2D Data first!!
 %This has side effect that it creates a figure if it is not yet present so needs to be after you create figure.
 % Swithches for user to play with ----------------------------------------
   SW_plotGR = 0;
   SW_plotTLS = 0;
   SW_f2plotlinTLS = 1;
 %-------------------------------------------------------------------------
 
 for kidn = kidn_iter % iterate over CKIDs.
    % IndexPsort is a vector of all the indices of powers for a specific
    % KID in Increasing order. 
    power_iter = IndexPopt(kidn); %[Index] - Popt -Number of P_read values. For now I look at 1 power. 
    %power_iter = IndexPsort{kidn,1}(7); %[Index] - All P_read Number of P_read values. For now I look at 1 power. 
    
    for p = power_iter % This is the power iterator over a specific kid 
    f(p) = figure('WindowState','maximized');
    ax(p) = axes('XScale','log','YScale','linear');
    hold(ax(p),'on')
    Tcolors = colormapJetJB(14); % this has been set to 14 such that you always generate the set of colors such that you always obtain the same color for a high temp.
        for nT = Tbath_iter
            
            % Fit TLS ~2 - 20Hz
            
            linDataY_TLS = NOISE(p).FFTnoise{nT}(:,4);
            TLS_coof_lin{p,nT} = LLS_TLS_SdB(freq(TLSf_min_i:TLSf_max_i),linDataY_TLS(TLSf_min_i:TLSf_max_i),Model_oneoverf,C_v_TLS_0);
            
            %f2 = figure;
            dBlogDataX = log10(freq(TLSf_min_i:TLSf_max_i));
            dBlogDataY = 10.*log10(NOISE(p).FFTnoise{nT}(:,4));
            %plot(dBlogDataX,dBlogDataY(TLSf_min_i:TLSf_max_i))
            TLS_coof_dBlog{p,nT} = polyfit(dBlogDataX,dBlogDataY(TLSf_min_i:TLSf_max_i),1);
            C_TLS = power(10,(TLS_coof_dBlog{p,nT}(2)/10));% This is the constant extracted from the linear fit in log space.
            gamma_TLS = -(TLS_coof_dBlog{p,nT}(1)/10); % "
            TLS_coof_dBlog_to_lin{p,nT} = [C_TLS gamma_TLS];
            figure(f(p)) % Brings the plotting figure back into focus.
            
            % Fit GR-noise + TLS ~400Hz - 100KHz
            
            linDataY_CB = NOISE(p).FFTnoise{nT}(:,4);
            linDataY_CB = linDataY_CB-Model_oneoverf([C_TLS gamma_TLS],freq); % subtract the TLS line..
            
            %Fitting GR noise spectrum 
            CB_coof_lin{p,nT} = LLS_CB_SdB(freq(CBf_min_i:CBf_max_i),linDataY_CB(CBf_min_i:CBf_max_i),Model_GR,C_v_CB_0,CB_lb,CB_ub);
            
            % Finding intersects
            fTLS_curr = @(x)abs(Model_oneoverf(TLS_coof_dBlog_to_lin{p,nT},x));
            fGR_curr = @(x)Model_GR(CB_coof_lin{p,nT},x);
            [f_inter(p,nT),S_inter(p,nT)] = findintersect_SdB(fTLS_curr,fGR_curr,fguess);
            plot(f_inter(p,nT),lintodb(S_inter(p,nT)),'o','Color',Tcolors(nT,:),'LineWidth',1.5,'HandleVisibility','off');
            %plot(1000,-170,'o','MarkerFaceColor', 'b','LineWidth',30)
            % Save these parameters in a struct/class?
            % contents: C_TLS, gamma, Res_TLS_gamma
            % C_GR, tau_qp, Residual_GR_Tau
            % P_read, Q_i, T_bath, Q_c , KID#
            % Measure for flatnes (Derivative!!)
            % Or you could say that a 2 db Diffence between GR-noise and measured data.
            
            % Fig. 1: Plotting the fractional frequency plot 
            %toplot = NOISE(p).FFTnoise{1}(:,4) > 0;
            %plot(NOISE(p).FFTnoise{1}(toplot,1),10*log10(NOISE(p).FFTnoise{1}(toplot,4)),...
            %         '-','color','k','LineWidth',2)
            toplot = NOISE(p).FFTnoise{nT}(:,4) > 0; % makes sure you only plot positive data. As you approach low noise you data can become negative due to noise.
            plot(NOISE(p).FFTnoise{nT}(toplot,1),10*log10(NOISE(p).FFTnoise{nT}(toplot,4)),...
                    '-','color',Tcolors(nT,:),'LineWidth',1.5)
                
            
            xlabel('F [Hz]');ylabel('S_F/F^2 [dBc/Hz]')
            %Old: %xlim([0.5,1e5]);grid on;ylim([-220,-140])
            xlim([0.1,5e5]);grid on;ylim([-220,-140])
            title(append('KID#',string(NOISE(p).KIDnumber)," |Power ",string(NOISE(p).ReadPower),"dBm"));
            %legendTvalues{nT+1-min(Tbath_iter)} = append(sprintf('%1.3f',NOISE(kidn).Temperature(nT)),' K');
            legendTvalues{nT} = append(sprintf('%1.3f',NOISE(kidn).Temperature(nT)),' K');
            
            %plot(freq(toplot),10*log10(Model_oneoverf(TLS_coof_lin{p,nT},freq(toplot))),'--','Color',Tcolors(nT,:),'LineWidth',1)
            if SW_plotTLS
            plot(freq(toplot),Model_line(TLS_coof_dBlog{p,nT},log10(freq(toplot))),'-.','Color','black','LineWidth',1,'HandleVisibility','off')
            end %Plotting TLS lines
            plot(freq(CBf_min_i:CBf_max_i),10.*log10(linDataY_CB(CBf_min_i:CBf_max_i)),'--','Color',Tcolors(nT,:),'LineWidth',0.5,'HandleVisibility','off');%
            %compensated lines!
            
            
            % GR lines
            if SW_plotGR
                plot(freq(toplot),10.*log10(Model_GR(CB_coof_lin{p,nT},freq(toplot))),'-.','Color','black','LineWidth',1,'HandleVisibility','off');
                %plot(freq(toplot),10.*log10(Model_GR(CB_coof_lin_CB{p,nT},freq(toplot))),'-.','Color',Tcolors(nT,:),'LineWidth',1,'HandleVisibility','off');
            end % Plots GR lines if user switch is high.
            toplotTLSplusGR =lintodb(dbtolin(Model_line(TLS_coof_dBlog{p,nT},log10(freq(toplot)))) + Model_GR(CB_coof_lin{p,nT},freq(toplot)));
            % Final model lines
            plot(freq(toplot),toplotTLSplusGR,'-.','Color',Tcolors(nT,:),'LineWidth',1.5,'HandleVisibility','off');
            xline(1/(2*pi*abs(CB_coof_lin{p,nT}(2))),'--','Color',Tcolors(nT,:),'LineWidth',0.25,'HandleVisibility','off')
        end % End Bath temperature.
    legendTvalues = legendTvalues(~cellfun('isempty',legendTvalues));
    % this removes the zero entries from the legendTvalues celll.
    hTl1(kidn) = legend(legendTvalues,'location','eastOutside');
    hTl1(kidn).ItemTokenSize = [10,10];
    xline(freq(TLSf_min_i),'--','Color','c','LineWidth',1,'HandleVisibility','off')
    xline(freq(TLSf_max_i),'--','Color','c','LineWidth',1,'HandleVisibility','off')
    xline(freq(CBf_min_i),'--','Color','m','LineWidth',1,'HandleVisibility','off')
    xline(freq(CBf_max_i),'--','Color','m','LineWidth',1,'HandleVisibility','off')
    
    hold(ax(p),'off')
    export_path_graph = append('../../../Export_Figures_noGit/LT254_DA_figures/Part_1_GR/f1KID',string(kidn),'Pread',sprintf('%2.0f',NOISE(p).ReadPower),'.png');
    exportgraphics(ax(p),export_path_graph)  
    end % End power 
end % End #KID

%% Part 1: Sec 2 | Plotting per Temperature + save figure!
% REMEMBER to load 2D Data first!!
 %This has side effect that it creates a figure if it is not yet present so needs to be after you create figure.
 % Swithches for user to play with ----------------------------------------
   SW_f2plotlinTLS = 0; % [Binary] plots the fit in linear space.
 %-------------------------------------------------------------------------
 for kidn = kidn_iter % iterate over CKIDs.
    % IndexPsort is a vector of all the indices of powers for a specific
    % KID in Increasing order. 
    power_iter = IndexPopt(kidn); %[Index] - Popt -Number of P_read values. For now I look at 1 power. 
    %power_iter = IndexPsort{kidn,1}(7); %[Index] - All P_read Number of P_read values. For now I look at 1 power. 
    
    for p = power_iter % This is the power iterator over a specific kid 
    
    Tcolors = colormapJetJB(14); % this has been set to 14 such that you always generate the set of colors such that you always obtain the same color for a high temp.
        for nT = Tbath_iter
            f2(nT) = figure('WindowState','maximized');
            ax2(nT) = axes('XScale','log','YScale','linear');
            hold(ax2(nT),'on')
            % Fit TLS ~2 - 20Hz
            %---TLS noise fit in linear space-----(Not used!)--------------
            linDataY_TLS = NOISE(p).FFTnoise{nT}(:,4);
            TLS_coof_lin{p,nT} = LLS_TLS_SdB(freq(TLSf_min_i:TLSf_max_i),linDataY_TLS(TLSf_min_i:TLSf_max_i),Model_oneoverf,C_v_TLS_0);
            if SW_f2plotlinTLS
                plot(freq(toplot),lintodb(Model_oneoverf(TLS_coof_lin{p,nT},freq(toplot))),'--','Color','black','LineWidth',1,'HandleVisibility','off');
                %plot(freq(toplot),10.*log10(Model_GR(CB_coof_lin_CB{p,nT},freq(toplot))),'-.','Color',Tcolors(nT,:),'LineWidth',1,'HandleVisibility','off');
            end
            %---TLS noise fit linear in Log-log space----------------------
            dBlogDataX = log10(freq(TLSf_min_i:TLSf_max_i));
            dBlogDataY = 10.*log10(NOISE(p).FFTnoise{nT}(:,4));
            
            TLS_coof_dBlog{p,nT} = polyfit(dBlogDataX,dBlogDataY(TLSf_min_i:TLSf_max_i),1);
            C_TLS = power(10,(TLS_coof_dBlog{p,nT}(2)/10));% This is the constant extracted from the linear fit in log space.
            gamma_TLS = -(TLS_coof_dBlog{p,nT}(1)/10); % "
           
            %-Fit GR-noise + TLS ~400Hz - 100KHz---------------------------
            
            linDataY_CB = NOISE(p).FFTnoise{nT}(:,4);
            linDataY_CB = linDataY_CB-Model_oneoverf([C_TLS gamma_TLS],freq); % subtract the TLS line..
            CB_coof_lin{p,nT} = LLS_CB_SdB(freq(CBf_min_i:CBf_max_i),linDataY_CB(CBf_min_i:CBf_max_i),Model_GR,C_v_CB_0,CB_lb,CB_ub);
            
            
            %---Finding intersects-----------------------------------------
            fTLS_curr = @(x)abs(Model_oneoverf(TLS_coof_dBlog_to_lin{p,nT},x));
            fGR_curr = @(x)Model_GR(CB_coof_lin{p,nT},x);
            [f_inter(p,nT),S_inter(p,nT)] = findintersect_SdB(fTLS_curr,fGR_curr,fguess);
            plot(f_inter(p,nT),lintodb(S_inter(p,nT)),'o','Color','black','LineWidth',1.5,'HandleVisibility','off');

            
            % Save these parameters in a struct/class?
            % contents: C_TLS, gamma, Res_TLS_gamma
            % C_GR, tau_qp, Residual_GR_Tau
            % P_read, Q_i, T_bath, Q_c , KID#
            % Measure for flatnes (Derivative!!)
            % Or you could say that a 2 db Diffence between GR-noise and measured data.
            
            % Fig. 1: Plotting the fractional frequency plot---------------
            toplot = NOISE(p).FFTnoise{nT}(:,4) > 0; % makes sure you only plot positive data. As you approach low noise you data can become negative due to noise.
            plot(NOISE(p).FFTnoise{nT}(toplot,1),10*log10(NOISE(p).FFTnoise{nT}(toplot,4)),...
                    '-','color',Tcolors(nT,:),'LineWidth',1.5)
                
            
            xlabel('F [Hz]');ylabel('S_F/F^2 [dBc/Hz]')
            %xlim([0.5,1e5]);grid on;ylim([-220,-140])
            xlim([0.1,5e5]);grid on;ylim([-220,-140])
            title(append('KID#',string(NOISE(p).KIDnumber)," |Power ",string(NOISE(p).ReadPower),"dBm"));
            legendTvalues2 = [{append(sprintf('%1.3f',NOISE(kidn).Temperature(nT)),' K')} {'GR model'} {'TLS model'} {'Combined model'} ];
            % TLS noise plot 
            plot(freq(toplot),Model_line(TLS_coof_dBlog{p,nT},log10(freq(toplot))),'-.','Color','black','LineWidth',1)
            plot(freq(CBf_min_i:CBf_max_i),10.*log10(linDataY_CB(CBf_min_i:CBf_max_i)),'--','Color',Tcolors(nT,:),'LineWidth',0.5,'HandleVisibility','off');%
            %compensated lines!
            
            
            % GR model plot
            plot(freq(toplot),10.*log10(Model_GR(CB_coof_lin{p,nT},freq(toplot))),':','Color','black','LineWidth',1.5);
            toplotTLSplusGR =lintodb(dbtolin(Model_line(TLS_coof_dBlog{p,nT},log10(freq(toplot)))) + Model_GR(CB_coof_lin{p,nT},freq(toplot)));
            % Total model: TLS + GR noise. 
            plot(freq(toplot),toplotTLSplusGR,'--','Color',Tcolors(nT,:),'LineWidth',1.5);
            hTl2(kidn) = legend(legendTvalues2,'location','eastOutside');
            %hTl2(kidn).ItemTokenSize = [20,10];
            xline(freq(TLSf_min_i),'--','Color','c','LineWidth',1,'HandleVisibility','off')
            xline(freq(TLSf_max_i),'--','Color','c','LineWidth',1,'HandleVisibility','off')
            xline(freq(CBf_min_i),'--','Color','m','LineWidth',1,'HandleVisibility','off')
            xline(freq(CBf_max_i),'--','Color','m','LineWidth',1,'HandleVisibility','off')
    
            hold(ax2(nT),'off')
            export_path_graph = append('../../../Export_Figures_noGit/LT254_DA_figures/Part_1_GR/f2KID',string(kidn),'T',sprintf('%1.3f',NOISE(p).Temperature(nT)),'.png');
            exportgraphics(ax2(nT),export_path_graph)
            %close(f2(nT)) % Closes the figures after use
        end % End Bath temperature. 
    end % End power 
end % End #KID


%% Part 2: TLS-noise coof.
% In this section we want to get plots of the noise at 1 KHz like in Gau
% paper.

% Loading the 120mK power data for maximum accuracy!! Gau paper is
% T=120mK
ChipInfo_path = ['..' filesep '..' filesep ]; %root path where data is, one higher than the scripts
FFT_power_subsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_120mK' filesep 'FFT' filesep 'Power'];% This is where the
load([ChipInfo_path FFT_power_subsubdir filesep 'NOISE_P.mat']) % 
freq = NOISE(1).FFTnoise{1,1}(:,1);





% Part 2: Sec 2
 % Swithches for user to play with ----------------------------------------
   SW_plotGR = 0;
   SW_plotTLS = 1;
   SW_f2plotlinTLS = 1;
 %-------------------------------------------------------------------------
 nT = 1;% choose index for temperature.
 for kidn = kidn_iter % iterate over CKIDs.
    % IndexPsort is a vector of all the indices of powers for a specific
    % KID in Increasing order. 
    %power_iter = IndexPopt(kidn); %[Index] - Popt -Number of P_read values. For now I look at 1 power. 
    power_iter = IndexPsort{kidn,1}; %[Index] - All P_read Number of P_read values. For now I look at 1 power. 
    
    for p = power_iter % This is the power iterator over a specific kid 
    f3(nT) = figure('WindowState','maximized');
    ax3(nT) = axes('XScale','log','YScale','linear');
    hold(ax3(nT),'on')
    Pcolors = colormapcoolSdB(max(power_iter)); % this has been set to 14 such that you always generate the set of colors such that you always obtain the same color for a high temp.
        
            
            % Fit TLS ~2 - 20Hz
            
            linDataY_TLS = NOISE(p).FFTnoise{nT}(:,4);
            TLS_coof_lin{p,nT} = LLS_TLS_SdB(freq(TLSf_min_i:TLSf_max_i),linDataY_TLS(TLSf_min_i:TLSf_max_i),Model_oneoverf,C_v_TLS_0);
            
            %f2 = figure;
            dBlogDataX = log10(freq(TLSf_min_i:TLSf_max_i));
            dBlogDataY = 10.*log10(NOISE(p).FFTnoise{nT}(:,4));
            %plot(dBlogDataX,dBlogDataY(TLSf_min_i:TLSf_max_i))
            TLS_coof_dBlog{p,nT} = polyfit(dBlogDataX,dBlogDataY(TLSf_min_i:TLSf_max_i),1);
            C_TLS = power(10,(TLS_coof_dBlog{p,nT}(2)/10));% This is the constant extracted from the linear fit in log space.
            gamma_TLS = -(TLS_coof_dBlog{p,nT}(1)/10); % "
            TLS_coof_dBlog_to_lin{p,nT} = [C_TLS gamma_TLS];
            figure(f3(nT)) % Brings the plotting figure back into focus.
            
            % Fit GR-noise + TLS ~400Hz - 100KHz
            
            linDataY_CB = NOISE(p).FFTnoise{nT}(:,4);
            linDataY_CB = linDataY_CB-Model_oneoverf([C_TLS gamma_TLS],freq); % subtract the TLS line..
            
            %Fitting GR noise spectrum 
            CB_coof_lin{p,nT} = LLS_CB_SdB(freq(CBf_min_i:CBf_max_i),linDataY_CB(CBf_min_i:CBf_max_i),Model_GR,C_v_CB_0,CB_lb,CB_ub);
            
            % Finding intersects
            fTLS_curr = @(x)abs(Model_oneoverf(TLS_coof_dBlog_to_lin{p,nT},x));
            fGR_curr = @(x)Model_GR(CB_coof_lin{p,nT},x);
            [f_inter(p,nT),S_inter(p,nT)] = findintersect_SdB(fTLS_curr,fGR_curr,fguess);
            %plot(f_inter(p,nT),lintodb(S_inter(p,nT)),'o','Color',Pcolors(p,:),'LineWidth',1.5,'HandleVisibility','off');
            %plot(1000,-170,'o','MarkerFaceColor', 'b','LineWidth',30)
            % Save these parameters in a struct/class?
            % contents: C_TLS, gamma, Res_TLS_gamma
            % C_GR, tau_qp, Residual_GR_Tau
            % P_read, Q_i, T_bath, Q_c , KID#
            % Measure for flatnes (Derivative!!)
            % Or you could say that a 2 db Diffence between GR-noise and measured data.
            
            % Fig. 1: Plotting the fractional frequency plot 
            %toplot = NOISE(p).FFTnoise{1}(:,4) > 0;
            %plot(NOISE(p).FFTnoise{1}(toplot,1),10*log10(NOISE(p).FFTnoise{1}(toplot,4)),...
            %         '-','color','k','LineWidth',2)
            toplot = NOISE(p).FFTnoise{nT}(:,4) > 0; % makes sure you only plot positive data. As you approach low noise you data can become negative due to noise.
            plot(NOISE(p).FFTnoise{nT}(toplot,1),10*log10(NOISE(p).FFTnoise{nT}(toplot,4)),...
                    '-','color',Pcolors(p,:),'LineWidth',1.5)
                
            
            xlabel('F [Hz]');ylabel('S_F/F^2 [dBc/Hz]')
            %Old: %xlim([0.5,1e5]);grid on;ylim([-220,-140])
            xlim([0.1,5e5]);grid on;ylim([-220,-140])
            title(append('KID#',string(NOISE(p).KIDnumber)," |T = ",string(NOISE(p).Temperature(nT)),"K"));
            %legendTvalues{nT+1-min(Tbath_iter)} = append(sprintf('%1.3f',NOISE(kidn).Temperature(nT)),' K');
            legendTvalues3{nT} = append(sprintf('%1.3f',NOISE(kidn).Temperature(nT)),' K');
            
            %plot(freq(toplot),10*log10(Model_oneoverf(TLS_coof_lin{p,nT},freq(toplot))),'--','Color',Pcolors(p,:),'LineWidth',1)
            if SW_plotTLS
            plot(freq(toplot),Model_line(TLS_coof_dBlog{p,nT},log10(freq(toplot))),'-.','Color','black','LineWidth',1,'HandleVisibility','off')
            end %Plotting TLS lines
            plot(freq(CBf_min_i:CBf_max_i),10.*log10(linDataY_CB(CBf_min_i:CBf_max_i)),'--','Color',Pcolors(p,:),'LineWidth',0.5,'HandleVisibility','off');%
            %compensated lines!
            
            
            % GR lines
            if SW_plotGR
                plot(freq(toplot),10.*log10(Model_GR(CB_coof_lin{p,nT},freq(toplot))),'-.','Color','black','LineWidth',1,'HandleVisibility','off');
                %plot(freq(toplot),10.*log10(Model_GR(CB_coof_lin_CB{p,nT},freq(toplot))),'-.','Color',Pcolors(p,:),'LineWidth',1,'HandleVisibility','off');
            end % Plots GR lines if user switch is high.
            toplotTLSplusGR =lintodb(dbtolin(Model_line(TLS_coof_dBlog{p,nT},log10(freq(toplot)))) + Model_GR(CB_coof_lin{p,nT},freq(toplot)));
            % Final model lines
            plot(freq(toplot),toplotTLSplusGR,'-.','Color',Pcolors(p,:),'LineWidth',1.5,'HandleVisibility','off');
            xline(1/(2*pi*abs(CB_coof_lin{p,nT}(2))),'--','Color',Pcolors(p,:),'LineWidth',0.25,'HandleVisibility','off')
    legendTvalues3 = legendTvalues3(~cellfun('isempty',legendTvalues3));
    % this removes the zero entries from the legendTvalues celll.
    hTl1(kidn) = legend(legendTvalues3,'location','eastOutside');
    hTl1(kidn).ItemTokenSize = [10,10];
    xline(freq(TLSf_min_i),'--','Color','c','LineWidth',1,'HandleVisibility','off')
    xline(freq(TLSf_max_i),'--','Color','c','LineWidth',1,'HandleVisibility','off')
    xline(freq(CBf_min_i),'--','Color','m','LineWidth',1,'HandleVisibility','off')
    xline(freq(CBf_max_i),'--','Color','m','LineWidth',1,'HandleVisibility','off')
    
    hold(ax3(nT),'off')
    export_path_graph = append('../../../Export_Figures_noGit/LT254_DA_figures/Part_1_GR/f3KID',string(kidn),'Pread',sprintf('%2.0f',NOISE(p).ReadPower),'.png');
    exportgraphics(ax3(nT),export_path_graph)  
    end % End power 
end % End #KID






























%% Appendix: (Old)TLS-noise coof.


% in this part of the script we want to do a TLS-noise study to find out
% how the TLS-noise depends on power. 
nT = 4 ;%Choose Temp index

 for kidn = kidn_iter % iterate over CKIDs.
    % IndexPsort is a vector of all the indices of powers for a specific
    % KID in Increasing order. 
    %power_iter = IndexPopt(kidn); %[Index] - Popt -Number of P_read values. For now I look at 1 power. 
    power_iter = IndexPsort{kidn,1}; %[Index] - All P_read Number of P_read values. For now I look at 1 power. 
    f3(kidn) = figure('WindowState','maximized');
    ax3(kidn) = axes('XScale','log','YScale','linear');
    hold(ax3(kidn),'on')
    for p = power_iter % This is the power iterator over a specific kid 
    
    Tcolors = colormapJetJB(14); % this has been set to 14 such that you always generate the set of colors such that you always obtain the same color for a high temp.
            
            % Fit TLS ~2 - 20Hz
            
            linDataY_TLS = NOISE(p).FFTnoise{nT}(:,4);
            TLS_coof_lin{p,nT} = LLS_TLS_SdB(freq(TLSf_min_i:TLSf_max_i),linDataY_TLS(TLSf_min_i:TLSf_max_i),Model_oneoverf,C_v_TLS_0);
            
            %f2 = figure;
            dBlogDataX = log10(freq(TLSf_min_i:TLSf_max_i));
            dBlogDataY = 10.*log10(NOISE(p).FFTnoise{nT}(:,4));
            %plot(dBlogDataX,dBlogDataY(TLSf_min_i:TLSf_max_i))
            TLS_coof_dBlog{p,nT} = polyfit(dBlogDataX,dBlogDataY(TLSf_min_i:TLSf_max_i),1);
            C_TLS = power(10,(TLS_coof_dBlog{p,nT}(2)/10));% This is the constant extracted from the linear fit in log space.
            gamma_TLS = -(TLS_coof_dBlog{p,nT}(1)/10); % "
            
            % Fit GR-noise + TLS ~400Hz - 100KHz
            
            linDataY_CB = NOISE(p).FFTnoise{nT}(:,4);
            linDataY_CB = linDataY_CB-Model_oneoverf([C_TLS gamma_TLS],freq); % subtract the TLS line..
            %f3 = figure
            %ax3 = axes('XScale','log','YScale','log');
            %plot(freq(CBf_min_i:CBf_max_i),linDataY_CB(CBf_min_i:CBf_max_i))
            %title('f3:Log-log')
            
            %C_TLS_corr = C_TLS./n_num
            CB_coof_lin{p,nT} = LLS_TLS_SdB(freq(CBf_min_i:CBf_max_i),linDataY_CB(CBf_min_i:CBf_max_i),Model_GR,C_v_CB_0);
            
            
            % Save these parameters in a struct/class?
            % contents: C_TLS, gamma, Res_TLS_gamma
            % C_GR, tau_qp, Residual_GR_Tau
            % P_read, Q_i, T_bath, Q_c , KID#
            % Measure for flatnes (Derivative!!)
            % Or you could say that a 2 db Diffence between GR-noise and measured data.
            
            % Fig. 1: Plotting the fractional frequency plot 
            %toplot = NOISE(p).FFTnoise{1}(:,4) > 0;
            %plot(NOISE(p).FFTnoise{1}(toplot,1),10*log10(NOISE(p).FFTnoise{1}(toplot,4)),...
            %         '-','color','k','LineWidth',2)
            toplot = NOISE(p).FFTnoise{nT}(:,4) > 0; % makes sure you only plot positive data. As you approach low noise you data can become negative due to noise.
            plot(NOISE(p).FFTnoise{nT}(toplot,1),10*log10(NOISE(p).FFTnoise{nT}(toplot,4)),...
                    '-','color',Tcolors(nT,:),'LineWidth',1.5)
                
            
            xlabel('F [Hz]');ylabel('S_F/F^2 [dBc/Hz]')
            %Old: %xlim([0.5,1e5]);grid on;ylim([-220,-140])
            xlim([0.1,5e5]);grid on;ylim([-220,-140])
            title(append('KID#',string(NOISE(p).KIDnumber)," |Power ",string(NOISE(p).ReadPower),"dBm"));
            %legendTvalues{nT+1-min(Tbath_iter)} = sprintf('%1.3f',NOISE(kidn).Temperature(nT));
            %plot(freq(toplot),10*log10(Model_oneoverf(TLS_coof_lin{p,nT},freq(toplot))),'--','Color',Tcolors(nT,:),'LineWidth',1)
            plot(freq(toplot),Model_line(TLS_coof_dBlog{p,nT},log10(freq(toplot))),'-.','Color',Tcolors(nT,:),'LineWidth',1,'HandleVisibility','off')
            plot(freq(CBf_min_i:CBf_max_i),10.*log10(linDataY_CB(CBf_min_i:CBf_max_i)),'--','Color',Tcolors(nT,:),'LineWidth',0.5,'HandleVisibility','off');%
            %compensated lines!
            
            
            
            plot(freq(toplot),10.*log10(Model_GR(CB_coof_lin{p,nT},freq(toplot))),'-.','Color',Tcolors(nT,:),'LineWidth',1,'HandleVisibility','off');
            toplotTLSplusGR =lintodb(dbtolin(Model_line(TLS_coof_dBlog{p,nT},log10(freq(toplot)))) + Model_GR(CB_coof_lin{p,nT},freq(toplot)));
            plot(freq(toplot),toplotTLSplusGR,'-.','Color',Tcolors(nT,:),'LineWidth',1.5,'HandleVisibility','off');
    
    %hTl1(kidn) = legend(legendTvalues,'location','eastOutside');
    %hTl1(kidn).ItemTokenSize = [10,10];
    xline(freq(TLSf_min_i),'--','Color','c','LineWidth',1,'HandleVisibility','off')
    xline(freq(TLSf_max_i),'--','Color','c','LineWidth',1,'HandleVisibility','off')
    xline(freq(CBf_min_i),'--','Color','m','LineWidth',1,'HandleVisibility','off')
    xline(freq(CBf_max_i),'--','Color','m','LineWidth',1,'HandleVisibility','off')
    
    hold(ax3(kidn),'off')
    export_path_graph = append('../../../Export_Figures_noGit/LT254_DA_figures/Part_1_GR/f1KID',string(kidn),'Pread',sprintf('%2.0f',NOISE(p).ReadPower),'.png');
    exportgraphics(ax3(kidn),export_path_graph)  
    end % End power 
end % End #KID












