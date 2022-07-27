% Last to run
% creates NEP plot (for now) at lowest T

clear variables;close all;clc

%================================================================================
% Input
%================================================================================
ChipInfo_path   = ['\\MARS\kid\KIDonSun\experiments\Entropy ADR\LT165_W2_Chip10' filesep]; %root path where data is, one higher than the scripts
FFTsubsubdir    = [ 'Noise_vs_T' filesep 'AllNoise' filesep '2D' ];
S21subsubdir    = [ 'S21' filesep 'Temp'];
S21file =       [S21subsubdir filesep 'ResponseS21.mat']; %full path
FFTfile =       [FFTsubsubdir filesep 'Noise_2D.mat'];
CsPSDfile   =   [FFTsubsubdir filesep 'CrossPSDNoise_2D'];
CsPSDfitfile=   [FFTsubsubdir filesep 'CrossPSDFit_2D','NOISEPdeplowT'];

tauinmax = 2e-3;
eta_pb = 0.4;
%================================================================================

addpath([pwd,filesep,'..',filesep,'subroutines']);

%load relevant data S21
load([ChipInfo_path,filesep,S21file],'KID');

%load relevant data noise files (level, lifetime etc)
load([ChipInfo_path,filesep,FFTfile],'NOISE','KIDnumbers','IndexP_sub_opt','IndexPopt' );
load([ChipInfo_path,filesep,CsPSDfile],'CrossPSDNOISE');
load([ChipInfo_path,filesep,CsPSDfitfile],'CrossPSDFit','NOISEPdeplowT');

ChipInfo.path = ChipInfo_path;clear ChipInfo_path; %needed to keep struct ok wrt path names when data is moved in between analysis ruyns

%define NEP folder
NEPfolder = [ChipInfo.path, filesep, 'NEP'] ;
if ~isfolder(NEPfolder)
    mkdir(NEPfolder);
end
%define output array
PoutData = zeros(length(KIDnumbers),11);

for kidn = 1:length(KIDnumbers) % LOOP OVER ALL UNIQUE KIDS,
    %construct filename
    legc = 1;%legend counter
    clear lstr

    kleur = colormapJetJB(length(NOISE(kidn).Temperature));
    kleur2 = colormapJetJB(length(IndexP_sub_opt{kidn}));
    
    for p=1:length(IndexP_sub_opt{kidn})% over Power
        KIDID = NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber;
        % response
        bla = [KID(:).KIDnumber];
        KIDindresp = (bla == KIDID);
        KIDtau_T = 1e-9*KID(KIDindresp).Ql./(pi*KID(KIDindresp).Fres); %tau res vs T in sec
        if sum(KIDindresp) ~= 1
            error('S21 files have multiple readout powers')
        end
        
        for nT=1:length(NOISE(IndexP_sub_opt{kidn}(p)).Temperature) % over T
            % create local variables with relevant data
            % Cross PSD results: NOte we also use the SR and Stheta
            % from the cross PSD program - the matlab is better than
            % labview
            F	 = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1); %frequency cross PSD
            CPSD = 10*log10(abs(-1*real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2))));      %cross PSD from V^2 to dBc/Hz
            SR = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,4);             %in dBc/Hz
            Stheta = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,3);         %in dBc/Hz
            %cross fit results
            tau = CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(nT);
            crosslevel = CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel(nT);
            crosssetup = CrossPSDFit(IndexP_sub_opt{kidn}(p)).setupnoise(nT);
            %P, T, Fres
            Pint = NOISE(IndexP_sub_opt{kidn}(p)).InternalPower(nT);
            Temperature = NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT);
            Tlegstr{nT} = ['T = ' num2str(1e3*Temperature,'%.0f')  ' mK'];
            Fres = NOISE(IndexP_sub_opt{kidn}(p)).Fres(nT);
           
            %get reponse and other params @ correct T by interpolation
            if min(KID(KIDindresp).Temperature) <= Temperature && max(KID(KIDindresp).Temperature) > Temperature
                dthetadN =  interp1(KID(KIDindresp).Temperature, KID(KIDindresp).ResponsivityM1(:,1), Temperature);
                dRdN =      interp1(KID(KIDindresp).Temperature, KID(KIDindresp).ResponsivityM1(:,2), Temperature);
                taures =    interp1(KID(KIDindresp).Temperature, KIDtau_T, Temperature);
            elseif min(KID(KIDindresp).Temperature) > Temperature && max(KID(KIDindresp).Temperature) > Temperature
                dthetadN =  KID(KIDindresp).ResponsivityM1(1,1);
                dRdN =      KID(KIDindresp).ResponsivityM1(1,2);
                taures =    KIDtau_T(1);
            else
                error('Repsonse T range terrible, cannot construct responses')
            end
            
            %store responsivity at lowest T for later. Note that this
            %parameter has no T dependence, but exact Tmin varies with
            %preadout possibly
            if Temperature == min(NOISE(IndexP_sub_opt{kidn}(p)).Temperature) %we are at min T
                NOISEPdeplowT(kidn).dthetadN_Tmin(p)    = dthetadN;
                NOISEPdeplowT(kidn).dRdN_Tmin(p)        = dRdN;
            end
            
            % Get the Noise levels devided by responsivity
            % CrossPSD: 1;F, 2=CPSD(V^2/Hz), 3=Stheta(dBc/Hz),
            % 4=SR(dBc/Hz), 5FITTed CPSD(V^2/Hz) renormlized
            % CrossPSD: 6=CPSD/(dthetadn*dRdN), 7=Stheta/dthetadn^2, 8=SR/dRdn^2,
            if CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel(nT) ~= 0 && ~isnan(CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel(nT))
                CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,5) = CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{nT}./(dthetadN.*dRdN);
            else
                CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,5) = zeros(size(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1)));
            end
            CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,6) = 10.^(CPSD/10)./(dthetadN.*dRdN); %i.e. V^2/Hz / resp^2
            CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,7) = 10.^(Stheta/10)./(dthetadN.*dthetadN);
            CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,8) = 10.^(SR/10)./(dRdN.*dRdN);
            NOISE(IndexP_sub_opt{kidn}(p)).dRdN(nT) = dRdN;
            NOISE(IndexP_sub_opt{kidn}(p)).dthetadN(nT) = dthetadN;
            crosslevel_norm = CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel(nT)./(dthetadN.*dRdN);%cross level is in V^2
            %crosssetup_norm = NOISE(IndexP_sub_opt{kidn}(p)).setupnoise(nT)./(dthetadN.*dRdN);
            CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel_norm(nT) = crosslevel_norm;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING Spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            figure(1e4*NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber + 1*abs(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower) );
            subplot(2,2,1) %Noise
            semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,4),...
                '-','color',kleur(nT,:),'LineWidth',1);%R
            hold on;grid on;xlim([1,1e5]);
            semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,3),...
                '-','color',kleur(nT,:),'LineWidth',2);%Theta
            
            subplot(2,2,2)
            loglog(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,6),...
                '-','color',kleur(nT,:),'LineWidth',2);hold on;grid on
            loglog(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,5),...
                '-k','LineWidth',1);%FIT
            
            subplot(2,2,3)
            loglog(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,7),...
                '-','color',kleur(nT,:),'LineWidth',2);hold on;grid on
            
            subplot(2,2,4)
            loglog(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,8),...
                '-','color',kleur(nT,:),'LineWidth',2);hold on;grid on
            
        end %end loop over all T
        
        %add legends and save
        subplot(2,2,1)
        legend('S_A','S_{\theta}');
        xlabel('F (Hz)');ylabel('S_x (dBc/Hz)')
        title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber) ]);
        
        subplot(2,2,2)
        xlabel('F [Hz]');ylabel('S_{Cross}/(dA/dN*d\theta/dN) [1/Hz]');xlim([0.5,1e5]);grid on;xlim([10,1e5]);
        title(['P_{read} = ',num2str(1*NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,'%.0f') ' dBm' ]);
            
        subplot(2,2,3)
        xlabel('F [Hz]');ylabel('S_\theta/(d\theta/dN)^2 [1/Hz]');xlim([0.5,1e5]);grid on;xlim([10,1e5]);
        
        subplot(2,2,4)
        xlabel('F [Hz]');ylabel('S_A/(dA/dN)^2 [1/Hz]');xlim([0.5,1e5]);grid on;xlim([10,1e5]);
        legend(Tlegstr);
        %save
        Figfile = [NEPfolder,filesep,...
            'KID',num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber,'%.0f'),'P_',num2str(-1*NOISE(IndexP_sub_opt{kidn}(p)).ReadPower,'%.0f'),'_crossPSD_norm_Tdep'];
        MakeGoodFigure(12,10,10,Figfile); %save png and fig
        close(gcf)
        
        %%%%%%%%%%%% NEP At low T for each power p %%%%%%%%%
        %% clear dangerous variables from T loop
        clear F	CPSD SR Stheta tau crosslevel crosssetup Pint Temperature Fres 
        [Tmin, Tminind] = min(NOISE(IndexP_sub_opt{kidn}(p)).Temperature);

        dthetadN    = NOISEPdeplowT(kidn).dthetadN_Tmin(p);
        dRdN        = NOISEPdeplowT(kidn).dRdN_Tmin(p);
        SR          = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,4);             %in dBc/Hz
        Stheta      = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,3);         %in dBc/Hz
        CPSD        = 10*log10(abs(-1*real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,2))));      %cross PSD from V^2 to dBc/Hz
        F           = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,1); %frequency cross PSD
        tau         = NOISEPdeplowT(kidn).Tau(p);
        taures      = interp1(KID(KIDindresp).Temperature, KIDtau_T, Tmin);
        crosslevel  = NOISEPdeplowT(kidn).crosslevel(p);
        %data from S21
        Tc          = KID(KIDindresp).Tc;           %not T dependent in the data!
        KIDID_resp  = KID(KIDindresp).KIDnumber;    %not T dependent in the data!
        Area        = KID(KIDindresp).Area;         %not T dependent in the data!
        Volume      = KID(KIDindresp).Volume;       %not T dependent in the data!
        Delta       = KID(KIDindresp).Delta;        %not T dependent in the data!
        Alu_L       = KID(KIDindresp).AluLength;    %not T dependent in the data!

        [NEPR,NEPtheta,NEPcross,Nqp,nqp,NEPGR,tau0] = ...
            getDARKNEP(dthetadN, dRdN, Stheta, SR, CPSD, F, eta_pb, tau, taures, crosslevel, Volume, Delta, Tc);
        NOISEPdeplowT(kidn).NEPR(p)
        NOISEPdeplowT(kidn).NEPtheta(p)
        NOISEPdeplowT(kidn).NEPcross(p)

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%% PLOT T dependence %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Final figure: T dependent sub-plots
        figure(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber );
        subplot(2,2,1) %tau vs T
        semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(:),CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau(:)*1e3,'--o',...
            'color',kleur2(p,:),'MarkerFaceColor',kleur2(p,:));hold on
        semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(:),CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumin(:)*1e3,'o','color',kleur2(p,:));
        semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(:),CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumax(:)*1e3,'o','color',kleur2(p,:));
        title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber) ]);
        
        subplot(2,2,2)%Noise 10 Hz normalized
        semilogy(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(:),CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel_norm(:),...
            '--o', 'color',kleur2(p,:),'MarkerFaceColor',kleur2(p,:));hold on
        legstr{p} = ['P_{Read} = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower)];
        MakeGoodFigure(12,10,10); %
        
        
    end %P
    
    subplot(2,2,1) %tau vs T
    xlabel('T (K)');ylabel('\tau (msec)');grid on;ylim([0.01 round(1.5*tauinmax*1000)]);
    
    subplot(2,2,2)%Noise 10 Hz normalized
    xlabel('T (K)');ylabel('C_{cross}/(d\theta/dN * dA/dN) @ 10 Hz');grid on;
    legend(legstr,'Location','SouthEast');
    title(['Alu length = ' num2str(Alu_L) 'nm']);
    
    %power dependent quantity at low T - renormalize v./(dthetadN.*dRdN)
    for p=1:length(IndexP_sub_opt{kidn})% over Power% over Power
        NOISEPdeplowT(kidn).crosslevel_norm(p) = NOISEPdeplowT(kidn).crosslevel(p)/( NOISEPdeplowT(kidn).dthetadN_Tmin(p)* NOISEPdeplowT(kidn).dRdN_Tmin(p));
    end
    %P dependent sub-plots
    subplot(2,2,3)
    semilogy(NOISEPdeplowT(kidn).ReadPower,NOISEPdeplowT(kidn).Tau*1e3,'--ro','MarkerFaceColor','r');hold on
    xlabel('P_{read} (dBm)');ylabel('\tau (msec)');grid on;%ylim([0.01 round(1.5*tauinmax*1000)]);xlim([0.1 0.35])
    title('\tau For all T < 145 mK')
    
    subplot(2,2,4)
    semilogy(NOISEPdeplowT(kidn).ReadPower,NOISEPdeplowT(kidn).crosslevel_norm,'--ro','MarkerFaceColor','r');hold on
    xlabel('P_{read} (dBm)');ylabel('S_{Cross} @ 10 Hz');grid on;%ylim([0.01 round(1.5*tauinmax*1000)]);xlim([0.1 0.35])
    title('C_{cross}/(d\theta/dN * dA/dN) For all T < 145 mK');
    MakeGoodFigure(14,12,10)
    
    Figfile2 = [NEPfolder,filesep,...
        'KID',num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber,'%.0f'),'_crossPSD_norm_Tdep_all_Pread'];
    MakeGoodFigure(12,10,10,Figfile2); %save png and fig
    
    
    
end %KID loop

save([NEPfolder,filesep,'NEP_2D.mat'],...
    'NOISE','KIDnumbers','IndexP_sub_opt','IndexPopt','CrossPSDFit','NOISEPdeplowT','CrossPSDNOISE');


rmpath([pwd,filesep,'..',filesep,'subroutines']);
