% NEP_GR_TdepV1
% reads td binary data, calulates and plots the PSD's
% saves the data into the relevant struct
%

clear variables;close all;clc

%================================================================================
% Input
%================================================================================
ChipInfo_path = [ '..' filesep '..'] ;%root path where data is, one higher than the scripts witgout filesep @ end;
S21file = [ 'S21' filesep 'Temp' filesep 'ResponseS21.mat'];
FFTfile = [ 'FFT' filesep 'Temp' filesep 'Noise_T.mat'];

eta_pb = 0.4;
lowT = 0.18;     %Tbelow tau is ~ constant.
FixTc = 1;  %1 = default. Tc not fitted in Kaplan fit. set 0 to fit it
%================================================================================

addpath([pwd,filesep,'..',filesep,'subroutines']);

%load relevant data S21
load([ChipInfo_path,filesep,S21file],'KID');

%load relevant data noise files (level, lifetime etc)
load([ChipInfo_path,filesep,FFTfile],'NOISE','CrossPSDNOISE','KIDnumbers');
ChipInfo.path = ChipInfo_path;clear ChipInfo_path; %needed to keep struct ok wrt path names when data is moved in between analysis ruyns

%define NEP folder
NEPfolder = [ChipInfo.path, filesep, 'NEP'] ;
if ~isfolder(NEPfolder)
    mkdir(NEPfolder);
end
%define output array
PoutData = zeros(length(KIDnumbers),11);

for kidn=1:length(KIDnumbers) % LOOP OVER ALL UNIQUE KIDS,
    %construct filename
    legc = 1;%legend counter
    clear lstr
    
    %indexing in response array and getting some parameters for easy
    %working
    KIDID = NOISE(kidn).KIDnumber;
    bla = [KID(:).KIDnumber]; 
    KIDindresp = (bla == KIDID);
    if sum(KIDindresp) == 1
        Fdesign = KID(KIDindresp).Fdesign;
        Delta = KID(KIDindresp).Delta;%not T dependent in the data!
        Tc = KID(KIDindresp).Tc;%not T dependent in the data!
        KIDID_resp = KID(KIDindresp).KIDnumber;%not T dependent in the data!
        Area = KID(KIDindresp).Area;%not T dependent in the data!
        Volume = KID(KIDindresp).Volume;%not T dependent in the data!
        Fres_resp = KID(KIDindresp).Fres(1);%Fres, not T dependent in the data!
        
        %==================================================================================================================================
        %PLOT
        %==================================================================================================================================
        figure(NOISE(kidn).KIDnumber)
        kleur = colormapJetJB(length(NOISE(kidn).Temperature));
        ki = 1;%index for color plots
        notauforthiskid = 1;
        for nT=1:length(NOISE(kidn).Temperature) % over T
            %check if we have a lifetome, othewsie we do not bother
            if NOISE(kidn).crosslevel(nT) ~= 0 %if this is 0. no tau fit, so no NEP
                % create local variables with relevant data
                % Cross PSD results: NOte we also use the SR and Stheta
                % from the cross PSD program - the matlab is better than
                % labview
                notauforthiskid = 0;
                F	 = CrossPSDNOISE(kidn).CrossPSD{nT}(:,1); %frequency cross PSD
                CPSD = 10*log10(abs(-1*real(CrossPSDNOISE(kidn).CrossPSD{nT}(:,2))));      %cross PSD in dBc/Hz
                SR = CrossPSDNOISE(kidn).CrossPSD{nT}(:,4);             %in dBc/Hz
                Stheta = CrossPSDNOISE(kidn).CrossPSD{nT}(:,3);         %in dBc/Hz
                %cross fit results
                tau = NOISE(kidn).tau(nT);
                crosslevel = NOISE(kidn).crosslevel(nT);
                crosssetup = NOISE(kidn).setupnoise(nT);
                %KID ids
                KIDID = NOISE(kidn).KIDnumber;
                Pint = NOISE(kidn).InternalPower(nT);
                Temperature = NOISE(kidn).Temperature(nT);
                Fres = NOISE(kidn).Fres(nT);
                % response
                KIDtau_T = 1e-9*KID(KIDindresp).Ql./(pi*KID(KIDindresp).Fres); %tau res vs T in sec
                
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
                    error('Repsonse T range terrible, cannot consruct NEP')
                end
                % Get the Noise levels devided by responsivity
                % CrossPSD: 1;F, 2=CPSD, 3=Stheta, 4=SR, 5FITTed CPSD
                % CrossPSD: 6=CPSD/(dthetadn*dRdN), 7=Stheta/dthetadn^2, 8=SR/dRdn^2, 
                % get the Noise levels re-normalised: convfert to V^2/Hz and devide by resp^2%

                CrossPSDNOISE(kidn).CrossPSD{nT}(:,6) = 10.^(CPSD/10)./(dthetadN.*dRdN); 
                CrossPSDNOISE(kidn).CrossPSD{nT}(:,7) = 10.^(Stheta/10)./(dthetadN.*dthetadN);
                CrossPSDNOISE(kidn).CrossPSD{nT}(:,8) = 10.^(SR/10)./(dRdN.*dRdN);
                NOISE(kidn).dRdN(nT) = dRdN;
                NOISE(kidn).dthetadN(nT) = dthetadN;

                
                % Get the NEP
                [NEPR,NEPtheta,NEPcross,Nqp,nqp,NEPGR,tau0] = getDARKNEP(dthetadN, dRdN, Stheta, SR, CPSD, F, eta_pb, tau, taures, crosslevel, Volume, Delta, Tc);
                
                %store NEP spectra
                NOISE(kidn).NEP{nT}(:,1) = CrossPSDNOISE(kidn).CrossPSD{nT}(:,1); % == F
                NOISE(kidn).NEP{nT}(:,2) = NEPR;       %Radius NEP
                NOISE(kidn).NEP{nT}(:,3) = NEPtheta;   %Phase NEP
                NOISE(kidn).NEP{nT}(:,4) = NEPcross;   %cross NEP
                NOISE(kidn).Nqp(nT)      = Nqp;
                NOISE(kidn).nqp(nT)      = nqp;
                NOISE(kidn).NEPGR(nT)    = NEPGR;
                NOISE(kidn).tau0_PSD(nT) = tau0;
                NOISE(kidn).dthetadN(nT) = dthetadN;
                NOISE(kidn).dRdN(nT)     = dRdN;
                
                % get NEP values smoothed above 10 Hz and below 1 kHz
                % NOISE(kidn).NEPRmin(nT,1) = F,
                % NOISE(kidn).NEPRmin(nT,2) = NEP value
                FR = find(NOISE(kidn).NEP{nT}(:,1) < 1e3 & NOISE(kidn).NEP{nT}(:,1) > 15);
                [NOISE(kidn).NEPRmin(nT,2),indi] = min(smooth(NEPR(FR),5));
                NOISE(kidn).NEPRmin(nT,1) = NOISE(kidn).NEP{nT}(FR(indi),1);
                [NOISE(kidn).NEPthetamin(nT,2),indi] = min(smooth(NEPtheta(FR),5));
                NOISE(kidn).NEPthetamin(nT,1) = NOISE(kidn).NEP{nT}(FR(indi),1);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                subplot(2,4,1) %Noise
                semilogx(CrossPSDNOISE(kidn).CrossPSD{nT}(:,1),CrossPSDNOISE(kidn).CrossPSD{nT}(:,4),...
                    '-','color',kleur(ki,:),'LineWidth',1);%R
                hold on;grid on;xlim([1,1e5]);
                semilogx(CrossPSDNOISE(kidn).CrossPSD{nT}(:,1),CrossPSDNOISE(kidn).CrossPSD{nT}(:,3),...
                    '-','color',kleur(ki,:),'LineWidth',2);%Theta
                
                %                 %Frequency Noise
                %                 subplot(2,4,2)
                %                 toplot = NOISE(kidn).FFTnoise{nT}(:,4) > 0;
                %                 semilogx(NOISE(kidn).FFTnoise{nT}(toplot,1),10*log10(NOISE(kidn).FFTnoise{nT}(toplot,4)),...
                %                     '-','color',kleur(ki,:),'LineWidth',1);hold on;grid on
                
                subplot(2,4,2)
                loglog(CrossPSDNOISE(kidn).CrossPSD{nT}(:,1),CrossPSDNOISE(kidn).CrossPSD{nT}(:,6),...
                    '-','color',kleur(ki,:),'LineWidth',1);hold on;grid on
                
                subplot(2,4,3)
                loglog(CrossPSDNOISE(kidn).CrossPSD{nT}(:,1),CrossPSDNOISE(kidn).CrossPSD{nT}(:,7),...
                    '-','color',kleur(ki,:),'LineWidth',1);hold on;grid on
                
                subplot(2,4,4)
                loglog(CrossPSDNOISE(kidn).CrossPSD{nT}(:,1),CrossPSDNOISE(kidn).CrossPSD{nT}(:,8),...
                    '-','color',kleur(ki,:),'LineWidth',1);hold on;grid on
                
                %                 subplot(2,4,3) %Cross PSD
                %                 semilogx(CrossPSDNOISE(kidn).CrossPSD{nT}(:,1),-1*real(CrossPSDNOISE(kidn).CrossPSD{nT}(:,2)),'-',...
                %                     'linewidth',2,'color',kleur(ki,:));hold on;grid on
                %                 semilogx(CrossPSDNOISE(kidn).CrossPSD{nT}(:,1),CrossPSDNOISE(kidn).CrossPSD{nT}(:,5),'-k');
                
                subplot(2,4,5)% NEP
                loglog(NOISE(kidn).NEP{nT}(:,1),NOISE(kidn).NEP{nT}(:,2),'-','Linewidth',1,'color',kleur(ki,:));%R
                grid on;hold on;
                loglog(NOISE(kidn).NEP{nT}(:,1),NOISE(kidn).NEP{nT}(:,3),'-','Linewidth',2,'color',kleur(ki,:));%theta
                loglog(NOISE(kidn).NEPRmin(nT,1),NOISE(kidn).NEPRmin(nT,2),'o','color',kleur(ki,:),'MarkerSize',8);
                loglog(NOISE(kidn).NEPthetamin(nT,1),NOISE(kidn).NEPthetamin(nT,2),'o','color',kleur(ki,:),'markerfacecolor',kleur(ki,:),'MarkerSize',8);
                
                
                subplot(2,4,6)%lifetime vs P
                semilogy(NOISE(kidn).Temperature(nT),NOISE(kidn).tau(nT)*1e3,'o',...
                    'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
                semilogy(NOISE(kidn).Temperature(nT),NOISE(kidn).taumin(nT)*1e3,'o','color',kleur(ki,:));
                semilogy(NOISE(kidn).Temperature(nT),NOISE(kidn).taumax(nT)*1e3,'o','color',kleur(ki,:));
                
                subplot(2,4,7)%nqp vs P Not made here
                
                subplot(2,4,8)%NEP vs P not made here
                
                legc = legc+2;%legend counter
                ki = ki +1; %kleur index
            else
                %store NEP values as NaN
                NOISE(kidn).Nqp(nT)      = NaN;
                NOISE(kidn).nqp(nT)      = NaN;
                NOISE(kidn).NEPGR(nT)    = NaN;
                NOISE(kidn).tau0_PSD(nT) = NaN;
                NOISE(kidn).dthetadN(nT) = NaN;
                NOISE(kidn).dRdN(nT)     = NaN;
            end %if statement if there is a fit
        end %end loop over all T
        
        %Now do T dependent stuff: Kaplan fit and Tau0 estimate
        
        if notauforthiskid == 0 % will nly be accessed for at least 1 fitted NEP
            %Kaplan fit
            if FixTc == 1
                [result,tau_fitted, nqpfullfit, nqp] = Fit_Kaplan2(NOISE(kidn).Temperature,NOISE(kidn).tau',lowT,Tc);%Tc
            elseif FixTc == 0
                [result,tau_fitted, nqpfullfit, nqp] = Fit_Kaplan2(NOISE(kidn).Temperature,NOISE(kidn).tau',lowT);%Tc
            else
                error('FixTc set wrong')
            end
            % I use nqp obtained from the Kaplan fit result tau_0 and measured lifetime.
            %rest plot
            NOISE(kidn).tau_0_Kaplan = result.tau_0;
            NOISE(kidn).NEPGR_Kaplan = 2*Delta/eta_pb * (Volume*nqp./(2*NOISE(kidn).tau')).^0.5;     %PdV 2.41; using the CORRECT single particle tau = 2*tau measured%
            NOISE(kidn).taufitted = tau_fitted;
            NOISE(kidn).nqp_Kaplan = nqp;
            NOISE(kidn).Tc_Kaplan = result.Tc;
            
            subplot(2,4,1)
            legend('S_A','S_{\theta}');
            xlabel('F (Hz)');ylabel('S_x (dBc/Hz)')
            title(['KID ' num2str(NOISE(kidn).KIDnumber) ]);
            
            subplot(2,4,2)
            xlabel('F [Hz]');ylabel('S_{Cross}/(dA/dN*d\theta/dN) [1/Hz]');xlim([0.5,1e5]);grid on;xlim([10,1e5]);
            
            subplot(2,4,3)
            xlabel('F [Hz]');ylabel('S_\theta/(d\theta/dN)^2 [1/Hz]');xlim([0.5,1e5]);grid on;xlim([10,1e5]);
                  
            subplot(2,4,4)
            xlabel('F [Hz]');ylabel('S_A/(dA/dN)^2 [1/Hz]');xlim([0.5,1e5]);grid on;xlim([10,1e5]);
            
            subplot(2,4,5)
            legend('NEP_A','NEP_{\theta}','Location','SouthEast');
            xlim([1,1e4]);ylim([1e-20 3e-17]);
            xlabel('F (Hz)');ylabel('NEP (W/\surd Hz)')
            title(['V = ' num2str(Volume) '\mum^3 , Area = ' num2str(Area) '\mum^2']);
            
            subplot(2,4,6)
            semilogy(NOISE(kidn).Temperature,1e3*NOISE(kidn).taufitted,'k-');
            xlabel('T  (K)');ylabel('\tau (msec)');grid on;ylim([0.01 3]);xlim([0.1 0.3])
            title(['Kaplan \tau_0 = ' num2str(1e9*NOISE(kidn).tau_0_Kaplan,'%.0f') ' nsec.']);
            
            subplot(2,4,7)
            semilogy(NOISE(kidn).Temperature,NOISE(kidn).nqp,'kx','MarkerSize',9);hold on;grid on;
            semilogy(NOISE(kidn).Temperature,NOISE(kidn).nqp_Kaplan,'ks','MarkerSize',9);%nqpfullfit
            semilogy(NOISE(kidn).Temperature,nqpfullfit,'k-');
            xlabel('T  (K)');ylabel('n_{qp} (\mum^{-3})');grid on;axis tight;
            xlim([0.1 0.3]);ylim([1 max(NOISE(kidn).nqp)*1.1]);
            legend(['n_{qp} from PSD_{cross}, Tc=' num2str(Tc ,'%0.2f') 'K'],...
                ['n_{qp} from \tau_0 from Kaplan fit and \tau_m. Tc=' num2str(NOISE(kidn).Tc_Kaplan ,'%0.2f') ' K'],...
                ['n_{qp} from \tau_0 and \tau from Kaplan fit. Tc=' num2str(NOISE(kidn).Tc_Kaplan ,'%0.2f') ' K']...
                ,'location','SouthEast')
            
            subplot(2,4,8)
            semilogy(NOISE(kidn).Temperature, NOISE(kidn).NEPGR,'kx','MarkerSize',9);hold on
            semilogy(NOISE(kidn).Temperature, NOISE(kidn).NEPGR_Kaplan,'ks','MarkerSize',9);
            ki = 1;
            for nT=1:length(NOISE(kidn).Temperature) % over T
                %check if we have a lifetome, othewsie we do not bother
                if NOISE(kidn).crosslevel(nT) ~= 0 %if this is 0. no tau fit, so no NEP
                    semilogy(NOISE(kidn).Temperature(nT),NOISE(kidn).NEPRmin(nT,2),'o',...
                        'color',kleur(ki,:));grid on;
                    semilogy(NOISE(kidn).Temperature(nT),NOISE(kidn).NEPthetamin(nT,2),'o',...
                        'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
                    ki = ki +1; %kleur index
                end
            end
            legend('NEP from from PSD_{cross}','NEP from FItted Kaplan \tau_0 +measured \tau',...
                'NEP_{R}','NEP_{\theta}','Location','NorthWest');
            xlabel('T  (K)');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-16]);xlim([0.1 0.3])
            title(['Fres = ' num2str(NOISE(kidn).Fres(1)) ' GHz']);
            
            
            Figfile = [NEPfolder,filesep,...
                'KID',num2str(NOISE(kidn).KIDnumber,'%.0f'),'_NEP_Tdep1'];
            MakeGoodFigure(16,15,9,Figfile); %save png and fig
        end % if tau is found from cross PSD
    else
        disp([num2str(KIDID) 'skipped'])
    end
end %KID loop

save([NEPfolder,filesep,'NEP_T.mat'],'NOISE','KIDnumbers');

rmpath([pwd,filesep,'..',filesep,'subroutines']);
