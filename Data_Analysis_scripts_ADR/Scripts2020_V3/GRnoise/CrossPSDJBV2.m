function CrossPSDJBV2
% reads td binary data, calulates and plots the PSD's
% saves the data into the relevant struct
%no further analysis
% NB:  % CrossPSD: 1;F, 2=CPSD(V^2/Hz), 3=Stheta(dBc/Hz),
                % 4=SR(dBc/Hz), 5FITTed CPSD(V^2/Hz)
%
clear variables;close all;clc

%================================================================================
% Input
ChipInfo_path = ['..' filesep '..' ]; %root path where data is, one higher than the scripts
PT2Ddep = 2;        %=0 for Tdep analysis, =1 for Pdep analysis, =2 for 2D analysis,=3 Only Popt(Sietse), error otheriwse
lowFupper = 500%300%200%70;    	%upper bound of low F range (level determined by quasiparticles) in Hz
lowFlower = 300%100%80%20;     %upper bound of low F range (level determined by quasiparticles)
nsigfitc = 1.2;   	%1.2 # sigma in the lowF noise that is condsidered enough to start a fit
plotall = 0;        %if =1 a plot of each fit is made and stored as png
closefigs = 0; % if 1 closes the P,T dep figures, leaves only the main ones (only 2D analysis)
tauinmax = 1.2e-3;    %maximum possible lifetime
highFpt = 2e4;    %1.3e4: high frequency reference pt for the setiup level for P dependence. ange will be 0.7 .. 1x this number%
LowTLIM     = 0.185; %all data below this is treated as lowT limit of lifetime and noise level
%================================================================================

if PT2Ddep == 0     %Tdep
    FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D_Popt'];                   %FFTsubdir = [filesep 'Noise_Powers_165mK' filesep 'FFT' filesep 'Power'];     %  
elseif PT2Ddep == 1 %Pdep
    FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D_Popt'];                   %FFTsubdir = [filesep 'Noise_Powers_165mK' filesep 'FFT' filesep 'Power'];     %
elseif PT2Ddep == 2 %2D
    %FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D_Popt' filesep 'CPSDMinusTLS'];
    FFTsubsubdir = ['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D'] ;  %FFTsubdir = [filesep 'Noise_Powers_165mK' filesep 'FFT' filesep 'Power'];     %
elseif PT2Ddep == 3 %Contour in T,P space across  Popt (Sietse)
    FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D_Popt' filesep 'CPSDMinusTLS'];                   %FFTsubdir = [filesep 'Noise_Powers_165mK' filesep 'FFT' filesep 'Power'];     %
else
    error('PT2Ddep not defined')
end



%================================================================================
if PT2Ddep == 0 %Tdep
    %Tdep
    matfile = 'Noise_T.mat';
    matfile2 = 'CrossPSDNoise_T';
    matfile3 = 'CrossPSDFit_T';
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile2],'CrossPSDNOISE');
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','KIDnumbers');
    CPSDfolder = [ChipInfo_path, filesep, 'CrossPSDTdep'] ;
elseif PT2Ddep == 1 %Pdep
    matfile = 'Noise_P.mat';
    matfile2 = 'CrossPSDNoise_P';
    matfile3 = 'CrossPSDFit_P';
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile2],'CrossPSDNOISE');
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','IndexP_sub_opt','KIDnumbers');
    CPSDfolder = [ChipInfo_path, filesep, 'CrossPSDPdep'] ;
elseif PT2Ddep == 2
    matfile = 'Noise_2D.mat';
    matfile2 = 'CrossPSDNoise_2D';
    matfile3 = 'CrossPSDFit_2D';
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','IndexP_sub_opt','KIDnumbers');
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile2],'CrossPSDNOISE');
    CPSDfolder = [ChipInfo_path, filesep, 'CrossPSD2D'] ;
elseif PT2Ddep == 3
    matfile = 'Noise_2D.mat';
    matfile2 = 'CrossPSDNoise_2D';
    matfile3 = 'CrossPSDFit_2D';
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','IndexP_sub_opt','KIDnumbers');
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile2],'CrossPSDNOISE');
    CPSDfolder = [ChipInfo_path, filesep, 'CrossPSD2D'] ;
else
    error('PT2Ddep not defined')
end


if ~isfolder(CPSDfolder)
    mkdir(CPSDfolder);
end

addpath([pwd,filesep,'..',filesep,'subroutines']);

ChipInfo.path = ChipInfo_path;clear ChipInfo_path; %needed to keep struct ok wrt path names when data is moved in between analysis ruyns
warning('off', 'MATLAB:Axes:NegativeDataInLogAxis')

for kidn= 1:length(KIDnumbers) % LOOP OVER ALL UNIQUE KIDS,
    %construct filename
    legc = 1;%legend counter
    clear lstr
    if PT2Ddep == 1
        kleur = colormapJetJB(length(IndexP_sub_opt{kidn}));
        for p=1:length(IndexP_sub_opt{kidn})% over Power
            %quick variables
            CPSD = abs(real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,2)));
            FCPSD = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,1);
            %auto condition to allow for fit
            lowF = mean(CPSD(lowFlower <= FCPSD & FCPSD <= lowFupper));
            s_lowF = std(CPSD(lowFlower < FCPSD & FCPSD <= lowFupper));
            highF = mean(CPSD(FCPSD > 0.7*highFpt & FCPSD < highFpt));
            %s_highF = std(CPSD(end-10:end));%ignored as being much smaller
            if plotall == 1
                %plot data for each
                figure(100000*NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber + abs(round(1000*NOISE(IndexP_sub_opt{kidn}(p)).ReadPower)));
                subplot(1,2,1)
                semilogx(FCPSD,CPSD,'b-');hold on;
                semilogx([lowFlower lowFupper],[lowF lowF],'r-','LineWidth',1);
                semilogx([lowFlower lowFupper],[lowF lowF]-s_lowF,'r--','LineWidth',1);semilogx([lowFlower lowFupper],[lowF lowF]+s_lowF,'r--','LineWidth',1);
                semilogx([1e4 2e4],[highF highF],'g-','LineWidth',2);
                subplot(1,2,2)
                loglog(FCPSD,CPSD,'b-');hold on;
                loglog([lowFlower lowFupper],[lowF lowF],'r-','LineWidth',1);
                loglog([lowFlower lowFupper],[lowF lowF]-s_lowF,'r--','LineWidth',1);semilogx([lowFlower lowFupper],[lowF lowF]+s_lowF,'r--','LineWidth',1);
                loglog([FCPSD(end-10) FCPSD(end)],[highF highF],'g-','LineWidth',2);
            end
            % see if we can fit and then do so
            if lowF - highF > nsigfitc*s_lowF
                % we fit starting at the low F range
                [tau,level,setupnoise,taumin,taumax] = crossfit(FCPSD(FCPSD > lowFlower),CPSD(FCPSD > lowFlower),lowF,tauinmax);
                CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{1} = level./(1 + (2*pi*tau*FCPSD).^2)+ setupnoise;
                if plotall == 1
                    subplot(1,2,1)
                    semilogx(FCPSD,CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{1},'k-');
                    legend('Data','GR level','GR-\sigma','GR+\sigma','highF level',['Fit, tau = ' num2str(tau*1000,'%0.3f') 'msec'])
                    subplot(1,2,2)
                    loglog(FCPSD,CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{1},'k-');
                end
            else
                tau = NaN;
                level = 0;
                setupnoise = 0;
                taumin = 0;taumax = 0;
            end
            
            %store results
            CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau = tau;
            CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumin = taumin;
            CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumax = taumax;
            CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel = level;
            CrossPSDFit(IndexP_sub_opt{kidn}(p)).setupnoise = setupnoise;
            if plotall == 1
                %finish plot
                subplot(1,2,1)
                xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
                grid on;axis tight;
                xlim([10,0.5e6]);
                subplot(1,2,2)
                xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
                grid on;axis tight;
                xlim([10,0.5e6]);
                title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber) ' @Pread = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower)...
                    ' dBm, T = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).Temperature)]);
                % save figure
                Figfile = [CPSDfolder,filesep,...
                    'KID',num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber,'%.0f'),'_',num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,'%.3g'),'dBm_',...
                    num2str(round(1000*NOISE(IndexP_sub_opt{kidn}(p)).Temperature)) 'mK_CrossPSD'];
                MakeGoodFigure(14,5,9,Figfile,1); %save png
                close gcf
            end
            
            %nice final figure
            figure(IndexP_sub_opt{kidn}(1))
            subplot(1,3,1) %Cross PSD linear
            semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,1),abs(real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,2))),'-',...
                'LineWidth',1,'color',kleur(p,:));hold on;
            if ~isnan(CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau)%add fit
                semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,1),CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{1},'-k');
                lstr{legc}=['P_{int} = ' num2str(round(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower)) ' dBm'];
                lstr{legc+1}=['\tau = ' num2str(round(1e6*CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau)/1e3) ' msec'];
                legc = legc+2;%legend counter
            else
                lstr{legc}=['P_{int} = ' num2str(round(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower)) ' dBm'];
                legc = legc+1;%legend counter
            end
            
            subplot(1,3,2) %Cross PSD log
            loglog(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,1),abs(real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,2))),'-',...
                'LineWidth',1,'color',kleur(p,:));hold on;
            if ~isnan(CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau)%add fit
                loglog(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,1),CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{1},'-k');
            end
            
            subplot(1,3,3)%lifetime vs P
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau*1e3,'o',...
                'color',kleur(p,:),'MarkerFaceColor',kleur(p,:));hold on
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumin*1e3,'o','color',kleur(p,:));
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumax*1e3,'o','color',kleur(p,:));
            
        end
        subplot(1,3,1)
        title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber) ' @ T = ' num2str(NOISE(IndexP_sub_opt{kidn}(1)).Temperature) ' K']);
        xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
        grid on;axis tight;
        xlim([10,0.5e6]);
        legend(lstr);
        
        subplot(1,3,2)
        xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
        grid on;axis tight;
        xlim([10,1e5]);
        
        subplot(1,3,3)
        xlabel('P_{Read} (dBm)');ylabel('\tau (msec)');grid on;ylim([0.01 2* tauinmax*1000]);
        Figfile = [CPSDfolder,filesep,...
            'KID',num2str(num2str(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber),'%.0f'),'_crossPSD_Pdep'];
        MakeGoodFigure(14,5,9,Figfile); %save png and fig
        save([ChipInfo.path,filesep,FFTsubsubdir,filesep,matfile3],'CrossPSDFit');
        
        %==================================================================================================================================
        %==================================================================================================================================
    elseif PT2Ddep == 0 %Tdepn(checked)
        %==================================================================================================================================
        %==================================================================================================================================
        kleur = colormapJetJB(length(NOISE(kidn).Temperature));
        
        for nT=1:length(NOISE(kidn).Temperature) % over T
            
            %quick variables
            CPSD = abs(real(CrossPSDNOISE(kidn).CrossPSD{nT}(:,2)));
            FCPSD = CrossPSDNOISE(kidn).CrossPSD{nT}(:,1);
            %auto condition to allow for fit
            lowF = mean(CPSD(lowFlower <= FCPSD & FCPSD <= lowFupper));
            s_lowF = std(CPSD(lowFlower < FCPSD & FCPSD <= lowFupper));
            highF = mean(CPSD(end-10:end));
            %s_highF = std(CPSD(end-10:end));%ignored as being much smaller
            if plotall == 1
                %plot data for each
                figure(100000*NOISE(kidn).KIDnumber + round(1000*NOISE(kidn).Temperature(nT)));
                subplot(1,2,1)
                semilogx(FCPSD,CPSD,'b-');hold on;
                semilogx([lowFlower lowFupper],[lowF lowF],'r-','LineWidth',1);
                semilogx([lowFlower lowFupper],[lowF lowF]-s_lowF,'r--','LineWidth',1);semilogx([lowFlower lowFupper],[lowF lowF]+s_lowF,'r--','LineWidth',1);
                semilogx([FCPSD(end-10) FCPSD(end)],[highF highF],'g-','LineWidth',2);
                subplot(1,2,2)
                loglog(FCPSD,CPSD,'b-');hold on;
                loglog([lowFlower lowFupper],[lowF lowF],'r-','LineWidth',1);
                loglog([lowFlower lowFupper],[lowF lowF]-s_lowF,'r--','LineWidth',1);semilogx([lowFlower lowFupper],[lowF lowF]+s_lowF,'r--','LineWidth',1);
                loglog([FCPSD(end-10) FCPSD(end)],[highF highF],'g-','LineWidth',2);
            end
            % see if we can fit and then do so
            if lowF - highF > nsigfitc*s_lowF
                % we fit starting at the low F range
                [tau,level,setupnoise,taumin,taumax] = crossfit(FCPSD(FCPSD > lowFlower),CPSD(FCPSD > lowFlower),lowF,tauinmax);
                CrossPSDFit(kidn).Fit{nT} = level./(1 + (2*pi*tau*FCPSD).^2)+ setupnoise;
                if plotall == 1
                    subplot(1,2,1)
                    semilogx(FCPSD,CrossPSDFit(kidn).Fit{nT},'k-');
                    legend('Data','GR level','GR-range','GR+range','highF level',['Fit, tau = ' num2str(tau*1000,'%0.3f') 'msec'])
                    subplot(1,2,2)
                    loglog(FCPSD,CrossPSDFit(kidn).Fit{nT},'k-');
                end
            else
                tau = NaN;
                level = 0;
                setupnoise = 0;
                taumin = 0;taumax = 0;
            end
            
            %store results
            CrossPSDFit(kidn).tau(nT) = tau;
            CrossPSDFit(kidn).taumin(nT) = taumin;
            CrossPSDFit(kidn).taumax(nT) = taumax;
            CrossPSDFit(kidn).crosslevel(nT) = level;
            CrossPSDFit(kidn).setupnoise(nT) = setupnoise;
            if plotall == 1
                %finish plot
                subplot(1,2,1)
                xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
                grid on;axis tight;
                xlim([10,0.5e6]);
                subplot(1,2,2)
                xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
                grid on;axis tight;
                xlim([10,0.5e6]);
                title(['KID ' num2str(NOISE(kidn).KIDnumber) ' @Pread = ' num2str(NOISE(kidn).ReadPower)...
                    ' dBm, T = ' num2str(NOISE(kidn).Temperature(nT))]);
                % save figure
                Figfile = [CPSDfolder,filesep,...
                    'KID',num2str(NOISE(kidn).KIDnumber,'%.0f'),'_',num2str(NOISE(kidn).ReadPower,'%.3g'),'dBm_',...
                    num2str(round(1000*NOISE(kidn).Temperature(nT))) 'mK_CrossPSD'];
                MakeGoodFigure(14,5,9,Figfile,1); %save png
                close gcf
            end
            
            %nice final figure
            figure(kidn)
            subplot(1,3,1)%noise + fit lin
            semilogx(CrossPSDNOISE(kidn).CrossPSD{nT}(:,1),abs(real(CrossPSDNOISE(kidn).CrossPSD{nT}(:,2))),'-',...
                'LineWidth',1,'color',kleur(nT,:));hold on;
            if ~isnan(CrossPSDFit(kidn).tau(nT))
                semilogx(CrossPSDNOISE(kidn).CrossPSD{nT}(:,1),CrossPSDFit(kidn).Fit{nT},'-k');
                lstr{legc}=['T = ' num2str(round(1000*NOISE(kidn).Temperature(nT))) 'mK'];
                lstr{legc+1}=['\tau = ' num2str(round(1e6*CrossPSDFit(kidn).tau(nT))/1e3) 'msec'];
                legc = legc + 2;
            else
                lstr{legc}=['T = ' num2str(round(1000*NOISE(kidn).Temperature(nT))) 'mK'];
                legc = legc + 1;
            end
            
            subplot(1,3,2)%noise + fit log
            loglog(CrossPSDNOISE(kidn).CrossPSD{nT}(:,1),abs(real(CrossPSDNOISE(kidn).CrossPSD{nT}(:,2))),'-',...
                'LineWidth',1,'color',kleur(nT,:));hold on;
            if ~isnan(CrossPSDFit(kidn).tau(nT))
                loglog(CrossPSDNOISE(kidn).CrossPSD{nT}(:,1),CrossPSDFit(kidn).Fit{nT},'-k');
            else
                lstr{legc}=['T = ' num2str(round(1000*NOISE(kidn).Temperature(nT))) 'mK'];
            end
            
            subplot(1,3,2)
            xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
            grid on;axis tight;
            xlim([10,1e5]);
            
            subplot(1,3,3)
            semilogy(NOISE(kidn).Temperature(nT),CrossPSDFit(kidn).tau(nT)*1e3,'o',...
                'color',kleur(nT,:),'MarkerFaceColor',kleur(nT,:));hold on
            semilogy(NOISE(kidn).Temperature(nT),CrossPSDFit(kidn).taumin(nT)*1e3,'o','color',kleur(nT,:));
            semilogy(NOISE(kidn).Temperature(nT),CrossPSDFit(kidn).taumax(nT)*1e3,'o','color',kleur(nT,:));
        end
        
        subplot(1,3,1)
        xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
        grid on;axis tight;
        xlim([10,0.5e6]);
        title(['KID ' num2str(NOISE(kidn).KIDnumber) ' @Pread = ' num2str(NOISE(kidn).ReadPower) ' dBm']);
        legend(lstr);
        
        subplot(1,3,3)
        xlabel('T (K)');ylabel('\tau (msec)');grid on;ylim([0.01 round(1.5*tauinmax*1000)]);xlim([0.1 0.35])
        Figfile = [CPSDfolder,filesep,...
            'KID',num2str(NOISE(kidn).KIDnumber,'%.0f'),'_crossPSD_Tdep'];
        MakeGoodFigure(15,5,9,Figfile); %save png and fig
        save([ChipInfo.path,filesep,FFTsubsubdir,filesep,matfile3],'CrossPSDFit');
        %==================================================================================================================================
        %==================================================================================================================================
    elseif PT2Ddep == 2 %2D(checked)
        %==================================================================================================================================
        %==================================================================================================================================
        kleur = colormapJetJB(length(NOISE(kidn).Temperature));
        kleur2 = colormapJetJB(length(IndexP_sub_opt{kidn}));
        for p=1:length(IndexP_sub_opt{kidn})% over Power
            for nT=1:length(NOISE(kidn).Temperature) % over T
                disp(['--> nT = ',string(nT)]);
                %quick variables
                CPSD = abs(real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2)));
                if length(CPSD) > 10%only do this if possible
                    FCPSD = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1);
                    %auto condition to allow for fit
                    lowF = mean(CPSD(lowFlower <= FCPSD & FCPSD <= lowFupper));
                    s_lowF = std(CPSD(lowFlower < FCPSD & FCPSD <= lowFupper));
                    highF = mean(CPSD(end-10:end));
                    %s_highF = std(CPSD(end-10:end));%ignored as being much smaller
                    if plotall == 1
                        %plot data for each
                        figure(1e7*NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber + 1e4*abs(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower) + round(1000*NOISE(kidn).Temperature(nT)));
                        subplot(1,2,1)
                        semilogx(FCPSD,CPSD,'b-');hold on;
                        semilogx([lowFlower lowFupper],[lowF lowF],'r-','LineWidth',1);
                        semilogx([lowFlower lowFupper],[lowF lowF]-s_lowF,'r--','LineWidth',1);semilogx([lowFlower lowFupper],[lowF lowF]+s_lowF,'r--','LineWidth',1);
                        semilogx([FCPSD(end-10) FCPSD(end)],[highF highF],'g-','LineWidth',2);
                        subplot(1,2,2)
                        loglog(FCPSD,CPSD,'b-');hold on;
                        loglog([lowFlower lowFupper],[lowF lowF],'r-','LineWidth',1);
                        loglog([lowFlower lowFupper],[lowF lowF]-s_lowF,'r--','LineWidth',1);semilogx([lowFlower lowFupper],[lowF lowF]+s_lowF,'r--','LineWidth',1);
                        loglog([FCPSD(end-10) FCPSD(end)],[highF highF],'g-','LineWidth',2);
                    end
                    % see if we can fit and then do so
                    if lowF - highF > nsigfitc*s_lowF
                        % we fit starting at the low F range
                        [tau,level,setupnoise,taumin,taumax] = crossfit(FCPSD(FCPSD > lowFlower),CPSD(FCPSD > lowFlower),lowF,tauinmax);
                        CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{nT} = level./(1 + (2*pi*tau*FCPSD).^2)+ setupnoise;
                        if plotall == 1
                            subplot(1,2,1)
                            semilogx(FCPSD,CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{nT},'k-');
                            legend('Data','GR level','GR-range','GR+range','highF level',['Fit, tau = ' num2str(tau*1000,'%0.3f') 'msec'])
                            subplot(1,2,2)
                            loglog(FCPSD,CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{nT},'k-');
                        end
                    else
                        tau = NaN;
                        level = 0;
                        setupnoise = 0;
                        taumin = 0;taumax = 0;
                    end
                    
                    %store results
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(nT) = tau;
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumin(nT) = taumin;
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumax(nT) = taumax;
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel(nT) = level;
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).setupnoise(nT) = setupnoise;
                    if plotall == 1
                        %finish plot
                        subplot(1,2,1)
                        xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
                        grid on;axis tight;
                        xlim([10,0.5e6]);
                        subplot(1,2,2)
                        xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
                        grid on;axis tight;
                        xlim([10,0.5e6]);
                        title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber) ' @Pread = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower)...
                            ' dBm, T = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT))]);
                        % save figure
                        Figfile = [CPSDfolder,filesep,...
                            'KID',num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber,'%.0f'),'_',num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,'%.3g'),'dBm_',...
                            num2str(round(1000*NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT))) 'mK_CrossPSD'];
                        MakeGoodFigure(14,5,9,Figfile,1); %save png
                        close gcf
                    end
                    
                    %nice final figure
                    
                    figure(1e4*NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber + 1*abs(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower) );
                    subplot(1,3,1)%noise + fit lin
                    semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),abs(real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2))),'-',...
                        'LineWidth',1,'color',kleur(nT,:));hold on;
                    if ~isnan(CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(nT))
                        semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{nT},'-k');
                        lstr{legc}=['T = ' num2str(round(1000*NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT))) 'mK'];
                        lstr{legc+1}=['\tau = ' num2str(round(1e6*CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(nT))/1e3) 'msec'];
                        legc = legc + 2;
                    else
                        lstr{legc}=['T = ' num2str(round(1000*NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT))) 'mK'];
                        legc = legc + 1;
                    end
                    
                    subplot(1,3,2)%noise + fit log
                    loglog(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),abs(real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2))),'-',...
                        'LineWidth',1,'color',kleur(nT,:));hold on;
                    if ~isnan(CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(nT))
                        loglog(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{nT},'-k');
                    else
                        lstr{legc}=['T = ' num2str(round(1000*NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT))) 'mK'];
                    end
                    
                    subplot(1,3,2)
                    xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
                    grid on;axis tight;
                    xlim([10,1e5]);
                    
                    subplot(1,3,3)
                    semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT),CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(nT)*1e3,'o',...
                        'color',kleur(nT,:),'MarkerFaceColor',kleur(nT,:));hold on
                    semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT),CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumin(nT)*1e3,'o','color',kleur(nT,:));
                    semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT),CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumax(nT)*1e3,'o','color',kleur(nT,:));
                    
                else
                    
                    %store results
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(nT) = NaN;
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumin(nT) = NaN;NOISE(IndexP_sub_opt{kidn}(p)).taumax(nT) = NaN;
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel(nT) = NaN;
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).setupnoise(nT) = NaN;
                end
            end
            %cont. figure per KID and per Pread
            warning('off','MATLAB:legend:IgnoringExtraEntries')
            subplot(1,3,1)
            xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
            grid on;axis tight;
            xlim([10,0.5e6]);
            title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber) ' @Pread = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower) ' dBm']);
            legend(lstr);
            
            subplot(1,3,3)
            xlabel('T (K)');ylabel('\tau (msec)');grid on;ylim([0.01 round(1.5*tauinmax*1000)]);xlim([0.1 0.35])
            Figfile = [CPSDfolder,filesep,...
                'KID',num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber,'%.0f'),'P_',num2str(-1*NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,'%.0f'),'_crossPSD_Tdep'];
            MakeGoodFigure(15,5,9,Figfile); %save png and fig
            if closefigs == 1
               close gcf 
            end
            %%%%%%%%%%%%%%%%%%
            
            
            %Final Fig
            figure(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber );
            subplot(1,3,1) %tau vs T
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(:),CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(:)*1e3,'--o',...
                'color',kleur2(p,:),'MarkerFaceColor',kleur2(p,:));hold on
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(:),CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumin(:)*1e3,'o','color',kleur2(p,:));
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(:),CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumax(:)*1e3,'o','color',kleur2(p,:));
            
            subplot(1,3,2)%Noise 10 Hz normalized 
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(:),CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel(:),...
                '--o', 'color',kleur2(p,:),'MarkerFaceColor',kleur2(p,:));hold on
            %legstr{p} = ['P_{Read} = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower)];
            legstr = ['P_{Read} = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower)];
            MakeGoodFigure(15,5,9); %
            
            subplot(1,3,3)%tau
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(1)*1e3,'--o',...
                'color',kleur2(p,:),'MarkerFaceColor',kleur2(p,:));hold on
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumin(1)*1e3,'o','color',kleur2(p,:));
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumax(1)*1e3,'o','color',kleur2(p,:));
           
            
        end
        subplot(1,3,1)
        xlabel('T (K)');ylabel('\tau (msec)');grid on;ylim([0.01 round(1.5*tauinmax*1000)]);xlim([0.1 0.35])
        title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber)]);
        
        subplot(1,3,2)
        xlabel('T (K)');ylabel('S_{Cross} @ 10 Hz');grid on;xlim([0.1 0.35])
        legend(legstr,'Location','Northwest');

        %create Pdependecie matrix. is NOT SAVED in the struct!
        for p=1:length(IndexP_sub_opt{kidn})% over Power% over Power
            Tlow = NOISE(IndexP_sub_opt{kidn}(1)).Temperature < LowTLIM;
            Tisnumb = ~isnan(CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau)';
            Ttouse = Tlow & Tisnumb;
            NOISEPdeplowT(kidn).ReadPower(p) = NOISE(IndexP_sub_opt{kidn}(p)).ReadPower;
            NOISEPdeplowT(kidn).Tau(p) = mean(CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(Ttouse));
            NOISEPdeplowT(kidn).crosslevel(p) = mean( CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel(Ttouse));
            NOISEPdeplowT(kidn).LowTLIM = LowTLIM;
        end
        
        subplot(1,3,3) %Pdep
        semilogy(NOISEPdeplowT(kidn).ReadPower,NOISEPdeplowT(kidn).Tau*1e3,'--ks','MarkerFaceColor','k');hold on
        xlabel('P_{Read} (dBm)');ylabel('\tau (msec)');grid on;%ylim([0.01 round(1.5*tauinmax*1000)]);xlim([0.1 0.35])
        title(['(Colors): \tau at lowest T. (Black): mean values below ' num2str(LowTLIM*1e3) ' mK '])
        Figfile2 = [CPSDfolder,filesep,...
            'KID',num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber,'%.0f'),'_CPSD_2D'];
        MakeGoodFigure(15,5,9,Figfile2); %save png and fig
        save([ChipInfo.path,filesep,FFTsubsubdir,filesep,matfile3],'CrossPSDFit');
    
        %==================================================================================================================================
        %==================================================================================================================================
    elseif PT2Ddep == 3 % Contour across Popt(checked)
        %==================================================================================================================================
        %==================================================================================================================================
        kleur = colormapJetJB(length(NOISE(kidn).Temperature));
        kleur2 = colormapJetJB(length(IndexP_sub_opt{kidn}));
        %for p=1:length(IndexP_sub_opt{kidn})% over Power
            for nT=1:length(NOISE(kidn).Temperature) % over T
                %Sietse Insert to find popt.
                % KiD: [1 2 3 4 5 6]
                P_Opt_order = [-94 -90 -94 -87 -92 -87];
                p = 1;
                while ~(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower == P_Opt_order(kidn))
                   p = p + 1; %What this does is it finds the right p for optimal power from the P_Opt_order list.
                end
                % Question how can i find index p such that
                % NOISE(IndexP_sub_opt{kidn}(p)).ReadPower is at P_opt
                % p should depend on kidn and nothing else
                %/Sietse
                %quick variables
                CPSD = abs(real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2)));
                if length(CPSD) > 10%only do this if possible
                    FCPSD = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1);
                    %auto condition to allow for fit
                    lowF = mean(CPSD(lowFlower <= FCPSD & FCPSD <= lowFupper));
                    s_lowF = std(CPSD(lowFlower < FCPSD & FCPSD <= lowFupper));
                    highF = mean(CPSD(end-10:end));
                    %s_highF = std(CPSD(end-10:end));%ignored as being much smaller
                    if plotall == 1
                        %plot data for each
                        figure(1e7*NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber + 1e4*abs(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower) + round(1000*NOISE(kidn).Temperature(nT)));
                        subplot(1,2,1)
                        semilogx(FCPSD,CPSD,'b-');hold on;
                        semilogx([lowFlower lowFupper],[lowF lowF],'r-','LineWidth',1);
                        semilogx([lowFlower lowFupper],[lowF lowF]-s_lowF,'r--','LineWidth',1);semilogx([lowFlower lowFupper],[lowF lowF]+s_lowF,'r--','LineWidth',1);
                        semilogx([FCPSD(end-10) FCPSD(end)],[highF highF],'g-','LineWidth',2);
                        subplot(1,2,2)
                        loglog(FCPSD,CPSD,'b-');hold on;
                        loglog([lowFlower lowFupper],[lowF lowF],'r-','LineWidth',1);
                        loglog([lowFlower lowFupper],[lowF lowF]-s_lowF,'r--','LineWidth',1);semilogx([lowFlower lowFupper],[lowF lowF]+s_lowF,'r--','LineWidth',1);
                        loglog([FCPSD(end-10) FCPSD(end)],[highF highF],'g-','LineWidth',2);
                    end
                    % see if we can fit and then do so
                    if lowF - highF > nsigfitc*s_lowF
                        % we fit starting at the low F range
                        [tau,level,setupnoise,taumin,taumax] = crossfit(FCPSD(FCPSD > lowFlower),CPSD(FCPSD > lowFlower),lowF,tauinmax);
                        CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{nT} = level./(1 + (2*pi*tau*FCPSD).^2)+ setupnoise;
                        if plotall == 1
                            subplot(1,2,1)
                            semilogx(FCPSD,CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{nT},'k-');
                            legend('Data','GR level','GR-range','GR+range','highF level',['Fit, tau = ' num2str(tau*1000,'%0.3f') 'msec'])
                            subplot(1,2,2)
                            loglog(FCPSD,CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{nT},'k-');
                        end
                    else
                        tau = NaN;
                        level = 0;
                        setupnoise = 0;
                        taumin = 0;taumax = 0;
                    end
                    
                    %store results
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(nT) = tau;
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumin(nT) = taumin;
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumax(nT) = taumax;
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel(nT) = level;
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).setupnoise(nT) = setupnoise;
                    if plotall == 1
                        %finish plot
                        subplot(1,2,1)
                        xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
                        grid on;axis tight;
                        xlim([10,0.5e6]);
                        subplot(1,2,2)
                        xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
                        grid on;axis tight;
                        xlim([10,0.5e6]);
                        title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber) ' @Pread = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower)...
                            ' dBm, T = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT))]);
                        % save figure
                        Figfile = [CPSDfolder,filesep,...
                            'KID',num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber,'%.0f'),'_',num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,'%.3g'),'dBm_',...
                            num2str(round(1000*NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT))) 'mK_CrossPSD'];
                        MakeGoodFigure(14,5,9,Figfile,1); %save png
                        close gcf
                    end
                    
                    %nice final figure
                    
                    figure(1e4*NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber + 1*abs(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower) );
                    subplot(1,3,1)%noise + fit lin
                    semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),abs(real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2))),'-',...
                        'LineWidth',1,'color',kleur(nT,:));hold on;
                    if ~isnan(CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(nT))
                        semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{nT},'-k');
                        lstr{legc}=['T = ' num2str(round(1000*NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT))) 'mK'];
                        lstr{legc+1}=['\tau = ' num2str(round(1e6*CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(nT))/1e3) 'msec'];
                        legc = legc + 2;
                    else
                        lstr{legc}=['T = ' num2str(round(1000*NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT))) 'mK'];
                        legc = legc + 1;
                    end
                    
                    subplot(1,3,2)%noise + fit log
                    loglog(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),abs(real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2))),'-',...
                        'LineWidth',1,'color',kleur(nT,:));hold on;
                    %disp(['--> nT = ',string(nT)]);
                    if ~isnan(CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(nT))
                        loglog(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{nT},'-k');
                    else
                        lstr{legc}=['T = ' num2str(round(1000*NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT))) 'mK'];
                    end
                    
                    subplot(1,3,2)
                    xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
                    grid on;axis tight;
                    xlim([10,1e5]);
                    
                    subplot(1,3,3)
                    semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT),CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(nT)*1e3,'o',...
                        'color',kleur(nT,:),'MarkerFaceColor',kleur(nT,:));hold on
                    semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT),CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumin(nT)*1e3,'o','color',kleur(nT,:));
                    semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT),CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumax(nT)*1e3,'o','color',kleur(nT,:));
                    
                else
                    
                    %store results
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(nT) = NaN;
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumin(nT) = NaN;NOISE(IndexP_sub_opt{kidn}(p)).taumax(nT) = NaN;
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel(nT) = NaN;
                    CrossPSDFit(IndexP_sub_opt{kidn}(p)).setupnoise(nT) = NaN;
                end
            end
            %cont. figure per KID and per Pread
            warning('off','MATLAB:legend:IgnoringExtraEntries')
            subplot(1,3,1)
            xlabel('F [Hz]');ylabel('S_{cross} [1/Hz]');
            grid on;axis tight;
            xlim([10,0.5e6]);
            title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber) ' @Pread = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower) ' dBm']);
            legend(lstr);
            
            subplot(1,3,3)
            xlabel('T (K)');ylabel('\tau (msec)');grid on;ylim([0.01 round(1.5*tauinmax*1000)]);xlim([0.1 0.35])
            Figfile = [CPSDfolder,filesep,...
                'KID',num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber,'%.0f'),'P_',num2str(-1*NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,'%.0f'),'_crossPSD_Tdep'];
            MakeGoodFigure(15,5,9,Figfile); %save png and fig
            if closefigs == 1
               close gcf 
            end
            %%%%%%%%%%%%%%%%%%
            
            
            %Final Fig
            figure(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber );
            subplot(1,3,1) %tau vs T
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(:),CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(:)*1e3,'--o',...
                'color',kleur2(p,:),'MarkerFaceColor',kleur2(p,:));hold on
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(:),CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumin(:)*1e3,'o','color',kleur2(p,:));
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(:),CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumax(:)*1e3,'o','color',kleur2(p,:));
            
            subplot(1,3,2)%Noise 10 Hz normalized 
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(:),CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel(:),...
                '--o', 'color',kleur2(p,:),'MarkerFaceColor',kleur2(p,:));hold on
            legstr{p} = ['P_{Read} = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower)];
            MakeGoodFigure(15,5,9); %
            
            subplot(1,3,3)%tau
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(1)*1e3,'--o',...
                'color',kleur2(p,:),'MarkerFaceColor',kleur2(p,:));hold on
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumin(1)*1e3,'o','color',kleur2(p,:));
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumax(1)*1e3,'o','color',kleur2(p,:));
           
            
        %end
        subplot(1,3,1)
        xlabel('T (K)');ylabel('\tau (msec)');grid on;ylim([0.01 round(1.5*tauinmax*1000)]);xlim([0.1 0.35])
        title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber)]);
        
        subplot(1,3,2)
        xlabel('T (K)');ylabel('S_{Cross} @ 10 Hz');grid on;xlim([0.1 0.35])
        legend(legstr{p},'Location','Northwest');

        %create Pdependecie matrix. is NOT SAVED in the struct!
        for p=1:length(IndexP_sub_opt{kidn})% over Power% over Power
            Tlow = NOISE(IndexP_sub_opt{kidn}(1)).Temperature < LowTLIM;
            Tisnumb = ~isnan(CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau)';
            Ttouse = Tlow & Tisnumb;
            NOISEPdeplowT(kidn).ReadPower(p) = NOISE(IndexP_sub_opt{kidn}(p)).ReadPower;
            NOISEPdeplowT(kidn).Tau(p) = mean(CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(Ttouse));
            NOISEPdeplowT(kidn).crosslevel(p) = mean( CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel(Ttouse));
            NOISEPdeplowT(kidn).LowTLIM = LowTLIM;
        end
        
        subplot(1,3,3) %Pdep
        semilogy(NOISEPdeplowT(kidn).ReadPower,NOISEPdeplowT(kidn).Tau*1e3,'--ks','MarkerFaceColor','k');hold on
        xlabel('P_{Read} (dBm)');ylabel('\tau (msec)');grid on;%ylim([0.01 round(1.5*tauinmax*1000)]);xlim([0.1 0.35])
        title(['(Colors): \tau at lowest T. (Black): mean values below ' num2str(LowTLIM*1e3) ' mK '])
        Figfile2 = [CPSDfolder,filesep,...
            'KID',num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber,'%.0f'),'_CPSD_2D'];
        MakeGoodFigure(15,5,9,Figfile2); %save png and fig
        save([ChipInfo.path,filesep,FFTsubsubdir,filesep,matfile3],'CrossPSDFit');
    end
    
    
end
disp(['Figures are in: ' Figfile])
rmpath([pwd,filesep,'..',filesep,'subroutines']);
