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
C_v_CB_0 = [(10^-18) 1.591549430918954e-06];% CTLS at 1 Hz , power with what  the power law decreases. TLS ~0.5

fguess = 2000; % [Hz] initial guess for the intersection point of the TLS and GR noise. 

% CB_lb = [10^-20 0.000000000000001];
% CB_ub = [10^-15 1];
CB_lb = [10^-20 (1/(2*pi*CBfitmax))];
CB_ub = [10^-15 (1/(2*pi*CBfitmin))];

%--Choosing your favorite KIDS,temperature and Reaout power(Moar Powa Baby)
kidn_iter = 2; %[Index] Select KID.
Tbath_iter = 8; %[Index] Number of Temperatures T_bath. Likely loop between 1:14;
SW_powerselector = 1 ;% choose power yourself or select Popt from file
power_selector = 7 ;% The index of the power that you want from low to high
%--End Choosing------------------------------------------------------------

%--Readout power selection ------------------------------------------------
if SW_powerselector % 1 case
f_power_iter = @(kidn)IndexPsort{kidn,1}(power_selector);
else % else
f_power_iter = @(kidn)IndexPopt(1,kidn);
end 
%-/Readout power selection -------------------------------------------------
tot_combi = length(concat_v(kidn_iter,IndexPsort,1));
f = gobjects(tot_combi,1);
ax = gobjects(tot_combi,1);
% Setup noise estimate
setup_noise_eval_freq = 350000;% [Hz]
[~,Snoise_i] = min(abs(freq-setup_noise_eval_freq)); % "

% Define fitting model 
% Here I will define the models that we will be using to fit the data. 
% This we will use to obtain a non-linear fit of the data. This gives heavy
% weight on the parts of the graph that are high on log log plot since here
% the difference is the biggest. 
Model_oneoverf = @(C_v,fdata)C_v(1).*power(fdata,-1*C_v(2)); %Model we use to fit C_v is the constants vector that we are trying to find.
Model_GR = @(C_v,fdata) C_v(1)./((1+power((2.*pi.*fdata.*C_v(2)),2).*(1))); %Model we use to fit C_v is the constants vector that we are trying to find.
Model_CB = @(C_v,fdata) C_v(1)./((1+power((2.*pi.*fdata.*C_v(2)),2).*(1)))+ Scurr_noise_level;
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
   SW_plotGR = 1;%[Binary] Plots the GR noise curve (Used)
   SW_plotTLS = 1;%[Binary] Plots polyfit in log log space (Used)
   SW_f2plotlinTLS = 1;%[Binary] plots the power law linear fit. Usually worse!
   SW_load2D = 1; %[Binary] Loads the data
   SW_plotsetup = 1;%[Binary] decides if you want to plot the found noise level
 %-------------------------------------------------------------------------
 if SW_load2D
load([ChipInfo_path FFTsubsubdir filesep 'NOISE_2D.mat'])
freq = NOISE(1).FFTnoise{1,1}(:,1);
end
 
 for kidn = kidn_iter % iterate over CKIDs.
    % IndexPsort is a vector of all the indices of powers for a specific
    % KID in Increasing order. 
    %power_iter = IndexPopt(kidn); %[Index] - Popt -Number of P_read values. For now I look at 1 power. 
    power_iter = f_power_iter(kidn); %[Index] - All P_read Number of P_read values. For now I look at 1 power. 
    
    for p = power_iter % This is the power iterator over a specific kid 
    
    f(p) = figure('WindowState','maximized');
    ax(p) = axes('XScale','log','YScale','linear');
    hold(ax(p),'on')
    Tcolors = colormapJetJB(14); % this has been set to 14 such that you always generate the set of colors such that you always obtain the same color for a high temp.
    
        for nT = Tbath_iter
            
            % Fit TLS ~2 - 20Hz
            
            linDataY_TLS = NOISE(p).FFTnoise{nT}(:,4);
            TLS_coof_lin{p,nT} = LLS_TLS_SdB(freq(TLSf_min_i:TLSf_max_i),linDataY_TLS(TLSf_min_i:TLSf_max_i),Model_oneoverf,C_v_TLS_0);
            S_noise_level(p,nT) = NOISE(p).FFTnoise{nT}(Snoise_i,4); % This is the setup noise level evaluated at 350000Hz.
            Scurr_noise_level = S_noise_level(p,nT) ;
            %Scurr_noise_level = dbtolin(-120)*(1/(4*NOISE(p).Ql(nT))^2);
            %Scurr_noise_level = 1E-25;
            % The anonymous function must be updated every time to
            % incorperate new function values
            Model_CB = @(C_v,fdata) C_v(1)./((1+power((2.*pi.*fdata.*C_v(2)),2).*(1)))+ Scurr_noise_level;

            if SW_plotsetup
                yline(lintodb(Scurr_noise_level))
            end
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
            CB_coof_lin{p,nT} = LLS_CB_SdB(freq(CBf_min_i:CBf_max_i),linDataY_CB(CBf_min_i:CBf_max_i),Model_CB,C_v_CB_0,CB_lb,CB_ub);
            
            % Finding intersects
            fTLS_curr = @(x)abs(Model_oneoverf(TLS_coof_dBlog_to_lin{p,nT},x));
            fGR_curr = @(x)Model_CB(CB_coof_lin{p,nT},x);
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
            %legendTvalues{nT+1-min(Tbath_iter)} = append(sprintf('%1.3f',NOISE(p).Temperature(nT)),' K');
            legendTvalues{nT} = append(sprintf('%1.3f',NOISE(p).Temperature(nT)),' K');
            
            %plot(freq(toplot),10*log10(Model_oneoverf(TLS_coof_lin{p,nT},freq(toplot))),'--','Color',Tcolors(nT,:),'LineWidth',1)
            if SW_plotTLS
            plot(freq(toplot),Model_line(TLS_coof_dBlog{p,nT},log10(freq(toplot))),'-.','Color','black','LineWidth',1,'HandleVisibility','off')
            end %Plotting TLS lines
            plot(freq(CBf_min_i:CBf_max_i),10.*log10(linDataY_CB(CBf_min_i:CBf_max_i)),'--','Color',Tcolors(nT,:),'LineWidth',0.5,'HandleVisibility','off');%
            %compensated lines!
            
            
            % GR lines
            if SW_plotGR
                plot(freq(toplot),10.*log10(Model_CB(CB_coof_lin{p,nT},freq(toplot))),'-.','Color','black','LineWidth',1,'HandleVisibility','off');
                %plot(freq(toplot),10.*log10(Model_CB(CB_coof_lin_CB{p,nT},freq(toplot))),'-.','Color',Tcolors(nT,:),'LineWidth',1,'HandleVisibility','off');
            end % Plots GR lines if user switch is high.
            toplotTLSplusGR =lintodb(dbtolin(Model_line(TLS_coof_dBlog{p,nT},log10(freq(toplot)))) + Model_CB(CB_coof_lin{p,nT},freq(toplot)));
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
disp('Part 1: Sec 1 - Done!')
%% Part 1: Sec 2 | Plotting per Temperature + save figure!
% REMEMBER to load 2D Data first!!
% Swithches for user to play with ----------------------------------------
   SW_f2plotlinTLS = 0; % [Binary] plots the fit in linear space.
   SW_load2D = 1; %[Binary] Loads the data
   SW_plotsetup = 1;%[Binary] decides if you want to plot the found noise level
 %-------------------------------------------------------------------------
if SW_load2D
load([ChipInfo_path FFTsubsubdir filesep 'NOISE_2D.mat'])
freq = NOISE(1).FFTnoise{1,1}(:,1);
end
 %This has side effect that it creates a figure if it is not yet present so needs to be after you create figure.
 
 for kidn = kidn_iter % iterate over CKIDs. 
    power_iter = f_power_iter(kidn); %[Index] - All P_read Number of P_read values. For now I look at 1 power. 
%     S_noise_level(p) = NOISE(p).FFTnoise{1}(Snoise_i,4); % This is the setup noise level evaluated at 350000Hz.
    for p = power_iter % This is the power iterator over a specific kid 
%     S_noise_level(p) = NOISE(p).FFTnoise{1}(Snoise_i,4); % This is the setup noise level evaluated at 350000Hz.
    Tcolors = colormapJetJB(14); % this has been set to 14 such that you always generate the set of colors such that you always obtain the same color for a high temp.
        for nT = Tbath_iter
            S_noise_level(p,nT) = NOISE(p).FFTnoise{nT}(Snoise_i,4); % This is the setup noise level evaluated at 350000Hz.
            Scurr_noise_level = S_noise_level(p,nT) ;
            %Scurr_noise_level = dbtolin(-120)*(1/(4*NOISE(p).Ql(nT))^2);
            %Scurr_noise_level = 1E-25;
            % The anonymous function must be updated every time to
            % incorperate new function values
            Model_CB = @(C_v,fdata) C_v(1)./((1+power((2.*pi.*fdata.*C_v(2)),2).*(1)))+ Scurr_noise_level;

            if SW_plotsetup
                yline(lintodb(Scurr_noise_level))
            end
            f2(nT) = figure('WindowState','maximized');
            ax2(nT) = axes('XScale','log','YScale','linear');
            hold(ax2(nT),'on')
%             Scurr_noise_level = S_noise_level(p) ;
            % Fit TLS ~2 - 20Hz
            %---TLS noise fit in linear space-----(Not used!)--------------
            linDataY_TLS = NOISE(p).FFTnoise{nT}(:,4);
            TLS_coof_lin{p,nT} = LLS_TLS_SdB(freq(TLSf_min_i:TLSf_max_i),linDataY_TLS(TLSf_min_i:TLSf_max_i),Model_oneoverf,C_v_TLS_0);
            if SW_f2plotlinTLS
                plot(freq(toplot),lintodb(Model_oneoverf(TLS_coof_lin{p,nT},freq(toplot))),'--','Color','black','LineWidth',1,'HandleVisibility','off');
                %plot(freq(toplot),10.*log10(Model_CB(CB_coof_lin_CB{p,nT},freq(toplot))),'-.','Color',Tcolors(nT,:),'LineWidth',1,'HandleVisibility','off');
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
            CB_coof_lin{p,nT} = LLS_CB_SdB(freq(CBf_min_i:CBf_max_i),linDataY_CB(CBf_min_i:CBf_max_i),Model_CB,C_v_CB_0,CB_lb,CB_ub);
            
            
            %---Finding intersects-----------------------------------------
            fTLS_curr = @(x)abs(Model_oneoverf(TLS_coof_dBlog_to_lin{p,nT},x));
            fGR_curr = @(x)Model_CB(CB_coof_lin{p,nT},x);
%             ftest(nT) = figure('WindowState','maximized');
%             axtest(nT) = axes('XScale','log','YScale','linear');
%             hold(axtest(nT),'on')
%             plot(freq,fTLS_curr(freq),freq,fGR_curr(freq) )
            [f_inter(p,nT),S_inter(p,nT)] = findintersect_SdB(fTLS_curr,fGR_curr,fguess);
            plot(f_inter(p,nT),lintodb(S_inter(p,nT)),'o','Color','black','LineWidth',1.5,'HandleVisibility','off');
%             hold(axtest(nT),'off')
            
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
            legendTvalues2 = [{append(sprintf('%1.3f',NOISE(p).Temperature(nT)),' K')} {'GR model'} {'TLS model'} {'Combined model'} ];
            % TLS noise plot 
            plot(freq(toplot),Model_line(TLS_coof_dBlog{p,nT},log10(freq(toplot))),'-.','Color','black','LineWidth',1)
            plot(freq(CBf_min_i:CBf_max_i),10.*log10(linDataY_CB(CBf_min_i:CBf_max_i)),'--','Color',Tcolors(nT,:),'LineWidth',0.5,'HandleVisibility','off');%
            %compensated lines!
            
            
            % GR model plot
            plot(freq(toplot),10.*log10(Model_CB(CB_coof_lin{p,nT},freq(toplot))),':','Color','black','LineWidth',1.5);
            toplotTLSplusGR =lintodb(dbtolin(Model_line(TLS_coof_dBlog{p,nT},log10(freq(toplot)))) + Model_CB(CB_coof_lin{p,nT},freq(toplot)));
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
clearvars legendPvalues;

SliceFreq = 1000;% [Hz] Freq at which we will take a slice of the spectrum and plot this according to Gau.
[~,Slice_i] = min(abs(freq-SliceFreq)); % "


% Part 2: Sec 2
 % Swithches for user to play with ----------------------------------------
   SW_plotGR = 0;
   SW_plotTLS = 1;
   SW_f2plotlinTLS = 1;
   SW_playJobs_done = 0; % Put zero to surpress Job's done sound
   SW_playGiorgio = 0; %Plays daft punk after your done. Because it is nice!
 %-------------------------------------------------------------------------
 nT = 1;% choose index for temperature.

 for kidn = kidn_iter % iterate over CKIDs.
    f3(kidn) = figure('WindowState','maximized');
    ax3(kidn) = axes('XScale','log','YScale','linear');
    hold(ax3(kidn),'on')
    % IndexPsort is a vector of all the indices of powers for a specific
    % KID in Increasing order. 
    %power_iter = IndexPopt(kidn); %[Index] - Popt -Number of P_read values. For now I look at 1 power. 
    power_iter = IndexPsort{kidn,1}; %[Index] - All P_read Number of P_read values. For now I look at 1 power. 
    
    for p = power_iter % This is the power iterator over a specific kid 
    
    %Pcolors = colormapcoolSdB(max(power_iter)); % this has been set to 14 such that you always generate the set of colors such that you always obtain the same color for a high temp.
    Pcolors = colormapcoolSdB(max(power_iter)-min(power_iter)+1)   
            
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
            figure(f3(kidn)) % Brings the plotting figure back into focus.
            
            % Fit GR-noise + TLS ~400Hz - 100KHz
            
            linDataY_CB = NOISE(p).FFTnoise{nT}(:,4);
            linDataY_CB = linDataY_CB-Model_oneoverf([C_TLS gamma_TLS],freq); % subtract the TLS line..
            
            %Fitting GR noise spectrum 
            CB_coof_lin{p,nT} = LLS_CB_SdB(freq(CBf_min_i:CBf_max_i),linDataY_CB(CBf_min_i:CBf_max_i),Model_CB,C_v_CB_0,CB_lb,CB_ub);
            
            % Finding intersects
            fTLS_curr = @(x)abs(Model_oneoverf(TLS_coof_dBlog_to_lin{p,nT},x));
            fGR_curr = @(x)Model_CB(CB_coof_lin{p,nT},x);
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
                    '-','color',Pcolors((p-min(power_iter)+1),:),'LineWidth',1.5)
            %Slice_i    
            plot(NOISE(p).FFTnoise{nT}(Slice_i,1),lintodb(NOISE(p).FFTnoise{nT}(Slice_i,4)),...
                    'o','color',Pcolors((p-min(power_iter)+1),:),'LineWidth',1,'MarkerFaceColor',Pcolors((p-min(power_iter)+1),:),'HandleVisibility','off')
            P2S2_SF_Eval_dB(kidn,p) = lintodb(NOISE(p).FFTnoise{nT}(Slice_i,4));
            P2S2_InternalPowerNaive(kidn,p) = NOISE(p).InternalPower;
            xlabel('F [Hz]');ylabel('S_F/F^2 [dBc/Hz]')
            %Old: %xlim([0.5,1e5]);grid on;ylim([-220,-140])
            xlim([0.1,5e5]);grid on;ylim([-220,-140])
            title(append('KID#',string(NOISE(p).KIDnumber)," | T = ",string(NOISE(p).Temperature(nT)),"K"));
            %legendTvalues{nT+1-min(Tbath_iter)} = append(sprintf('%1.3f',NOISE(p).Temperature(nT)),' K');
            legendPvalues{p} = append(sprintf('%3.0f',NOISE(p).ReadPower),' dBc');
            %plot(freq(toplot),10*log10(Model_oneoverf(TLS_coof_lin{p,nT},freq(toplot))),'--','Color',Pcolors(p,:),'LineWidth',1)
            if SW_plotTLS
            plot(freq(toplot),Model_line(TLS_coof_dBlog{p,nT},log10(freq(toplot))),'-.','Color','black','LineWidth',1,'HandleVisibility','off')
            end %Plotting TLS lines
            plot(freq(CBf_min_i:CBf_max_i),10.*log10(linDataY_CB(CBf_min_i:CBf_max_i)),'--','Color',Pcolors((p-min(power_iter)+1),:),'LineWidth',0.5,'HandleVisibility','off');%
            %compensated lines!
            
            
            % GR lines
            if SW_plotGR
                plot(freq(toplot),10.*log10(Model_CB(CB_coof_lin{p,nT},freq(toplot))),'-.','Color','black','LineWidth',1,'HandleVisibility','off');
                %plot(freq(toplot),10.*log10(Model_CB(CB_coof_lin_CB{p,nT},freq(toplot))),'-.','Color',Pcolors(p,:),'LineWidth',1,'HandleVisibility','off');
            end % Plots GR lines if user switch is high.
            toplotTLSplusGR =lintodb(dbtolin(Model_line(TLS_coof_dBlog{p,nT},log10(freq(toplot)))) + Model_CB(CB_coof_lin{p,nT},freq(toplot)));
            % Final model lines
            plot(freq(toplot),toplotTLSplusGR,'-.','Color',Pcolors((p-min(power_iter)+1),:),'LineWidth',1.5,'HandleVisibility','off');
            xline(1/(2*pi*abs(CB_coof_lin{p,nT}(2))),'--','Color',Pcolors((p-min(power_iter)+1),:),'LineWidth',0.25,'HandleVisibility','off')
    
    xline(freq(TLSf_min_i),'--','Color','c','LineWidth',1,'HandleVisibility','off')
    xline(freq(TLSf_max_i),'--','Color','c','LineWidth',1,'HandleVisibility','off')
    xline(freq(CBf_min_i),'--','Color','m','LineWidth',1,'HandleVisibility','off')
    xline(freq(CBf_max_i),'--','Color','m','LineWidth',1,'HandleVisibility','off')
    
    
     
    end % End power 
    

    legendPvalues = legendPvalues(~cellfun('isempty',legendPvalues));
    % this removes the zero entries from the legendTvalues celll.
    hTl1(kidn) = legend(legendPvalues,'location','eastOutside');
    hTl1(kidn).ItemTokenSize = [10,10];
    hold(ax3(kidn),'off')   
    export_path_graph = append('../../../Export_Figures_noGit/LT254_DA_figures/Part_2_TLS/f3KID',string(kidn),'Pread',sprintf('%2.0f',NOISE(p).ReadPower),'.png');
    exportgraphics(ax3(kidn),export_path_graph)
 
end % End #KID


%% Part 2. Sec 2 - Plotting the curves like in Gau
f4 = figure('WindowState','maximized');
ax4 = axes('XScale','linear','YScale','linear');
hold(ax4,'on')
% Info: NOISE(p).InternalPower = 10*log10((2/pi)*10.^(NOISE(p).ReadPower/10).*(NOISE(p).Ql.^2./NOISE(p).Qc));
% WARNING: This definition might not be correct. And according to Akira
% this can be a tricky topic to go into!

Customcolormap = [0 0 1  ; 0 0 0.9; 0 0 0.8 ;0 1 0 ;0 0.9 0;0 0.8 0];



for kidn=kidn_iter
p1 = plot(P2S2_InternalPowerNaive(kidn,IndexPsort{kidn,1}),P2S2_SF_Eval_dB(kidn,IndexPsort{kidn,1}),'-o');
set(p1,'Color',Customcolormap(kidn,:));
set(p1,'MarkerFaceColor',Customcolormap(kidn,:));
% This needs to be sorted in Power
% I need to add the lines from Gau paper

end
plot(linspace(-90,-30,2),Model_line([-0.5 -192],linspace(-80,-40,2)),'--','LineWidth',2,'Color','black')
f4legendstr = [{'\#1 ( C7G3: 2-2-2 ?)'} {'\#2 ( C9G3: 2-2-2 ?)'} {'\#3 ( C8G3: 2-2-2 ?)'} {'\#4 ( C10G4: 4-4-4 ?)'} {'\#5 ( C11G4: 4-4-4 ?)'} {'\#6 ( C12G4: 4-4-4 ?)'} {'$P_{int}^{-\frac{1}{2}}$ (Theory)'}];
f4lgd = legend(f4legendstr, 'interpreter','latex');
f4lgd.FontSize = 14;
hold(ax4,'off')
grid on
xlabel('P_{int}^{*}');ylabel('S_F/F^2 @1KHz [dBc/Hz]')

xlim([-110,-15]);grid on;ylim([-200,-150])
title('Gau compare - Need to add Gau data!')

export_path_graph = append('../../../Export_Figures_noGit/LT254_DA_figures/Part_2_TLS/f4Gau','.png');
exportgraphics(ax4,export_path_graph)

disp("Jobs done!!")
if SW_playJobs_done
[y, Fs] = audioread('Jobs_Done.mp3');
player = audioplayer(y, Fs);
play(player);
end

if SW_playGiorgio 
[y, Fs] = audioread('Giorgio_by_moroder.mp3');
player = audioplayer(y, Fs);
play(player);
disp('My name is Giovanni Giorgio. But everybody calls me Giorgio!')
disp('Click any button to stop music!')
pause
stop(player)
end

%% Appendix: (Old)TLS-noise coof.





%% Test 1: So if we make an error in the gamma coof. does this explain the discrapancy?















