% NEP_GR_PdepV1
% reads td binary data, calulates and plots the PSD's
% saves the data into the relevant struct
%

clear variables;close all;clc

%================================================================================
% Input
%================================================================================
ChipInfo_path   = ['\\MARS\kid\KIDonSun\experiments\Entropy ADR\LT197-chip9' filesep]; %root path where data is, one higher than the scripts
FFTsubsubdir    = [ 'Noise 120mK' filesep 'FFT' filesep 'Power' ];
S21subsubdir    = [ 'S21' filesep 'Temp'];
S21file =       [S21subsubdir filesep 'ResponseS21.mat']; %full path
FFTfile =       [FFTsubsubdir filesep 'Noise_P.mat'];
CsPSDfile   =   [FFTsubsubdir filesep 'CrossPSDNoise_P'];
CsPSDfitfile=   [FFTsubsubdir filesep 'CrossPSDFit_P'];

eta_pb = 0.5;
%tau_0 = 1500e-9; %in sec, estimated or obtained from T dep data. Could be KID:KID coded (data is present) but getting good data is sketchy and this parameter shold be the same for all devices.
%================================================================================

addpath([pwd,filesep,'..',filesep,'subroutines']);

%load relevant data S21
load([ChipInfo_path,filesep,S21file],'KID');

%load relevant data noise files (level, lifetime etc)
load([ChipInfo_path,filesep,FFTfile],'NOISE','KIDnumbers','IndexP_sub_opt','IndexPopt','IndexPsort');
load([ChipInfo_path,filesep,CsPSDfile],'CrossPSDNOISE');
load([ChipInfo_path,filesep,CsPSDfitfile],'CrossPSDFit');

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
    ki = 1;%index for color plots
    notauforthiskid = 1;
    for p=1:length(IndexP_sub_opt{kidn})% over Power
        KIDID = NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber;
        % response
        bla = [KID(:).KIDnumber];
        KIDindresp = (bla == KIDID);
        if sum(KIDindresp) ~= 1
            error('S21 files have multiple readout powers')
        end
        Fdesign = KID(KIDindresp).Fdesign;
        Delta = KID(KIDindresp).Delta;%not T dependent in the data!
        Tc = KID(KIDindresp).Tc;%not T dependent in the data!
        KIDID_resp = KID(KIDindresp).KIDnumber;%not T dependent in the data!
        Area = KID(KIDindresp).Area;%not T dependent in the data!
        Volume = KID(KIDindresp).Volume;%not T dependent in the data!
        Fres_resp = KID(KIDindresp).Fres(1);%Fres, not T dependent in the data!
        AluLength = KID(KIDindresp).AluLength;
        AluV = KID(KIDindresp).Volume;
        dxdN = KID(KIDindresp).ddxdNqp(1);

        figure(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber)
        kleur = colormapJetJB(sum([CrossPSDFit(IndexP_sub_opt{kidn}).crosslevel]' ~= 0));
        
        %check if we have a lifetome, othewsie we do not bother
        if CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel ~= 0 %not 0, so we have tua fit and can continue
            % create local variables with relevant data
            % Cross PSD results: NOte we also use the SR and Stheta
            % from the cross PSD program - the matlab is better than
            % labview
            notauforthiskid = 0;
            F	 = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,1); %frequency cross PSD
            CPSD = 10*log10(abs(-1*real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,2)))); %cross PSD
            SR = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,4);
            Stheta = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,3);
            %cross fit results
            tau = CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau;
            crosslevel = CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel;
            crosssetup = CrossPSDFit(IndexP_sub_opt{kidn}(p)).setupnoise;
            %KID ids
            KIDIDcheck = NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber;
            if KIDID ~= KIDIDcheck
                error('KID ID problems!!')
            end
            Pint = NOISE(IndexP_sub_opt{kidn}(p)).InternalPower;
            Temperature = NOISE(IndexP_sub_opt{kidn}(p)).Temperature;
            Fres = NOISE(IndexP_sub_opt{kidn}(p)).Fres;
            % response
            KIDtau_T = 1e-9*KID(KIDindresp).Ql./(pi*KID(KIDindresp).Fres); %tau res vs T in sec
            
            %disp(['Check Fres: ' num2str([Fres Fres_resp])]);
            %disp(['Check ID: ' num2str([KIDID KIDID_resp])]);
            
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
            % get the Noise levels re-normalised: convfert to V^2/Hz and
            % devide by resp^2% 
            CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,5) = CrossPSDFit(IndexP_sub_opt{kidn}(p)).Fit{1}./(dthetadN.*dRdN);
            CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,6) = 10.^(CPSD/10)./(dthetadN.*dRdN);
            CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,7) = 10.^(Stheta/10)./(dthetadN.*dthetadN);
            CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,8) = 10.^(SR/10)./(dRdN.*dRdN);
            
            NOISE(IndexP_sub_opt{kidn}(p)).dRdN = dRdN;
            NOISE(IndexP_sub_opt{kidn}(p)).dthetadN = dthetadN;
            
            % Get the NEP
            [NEPR,NEPtheta,NEPcross,Nqp,nqp,NEPGR,tau0] = getDARKNEP(dthetadN, dRdN, Stheta, SR, CPSD, F, eta_pb, tau, taures, crosslevel, Volume, Delta, Tc);
            
            NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,1) = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,1); % == F
            NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,2) = NEPR;
            NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,3) = NEPtheta;
            NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,4) = NEPcross;
            NOISE(IndexP_sub_opt{kidn}(p)).Nqp = Nqp;
            NOISE(IndexP_sub_opt{kidn}(p)).nqp = nqp;
            NOISE(IndexP_sub_opt{kidn}(p)).NEPGR = NEPGR;
            NOISE(IndexP_sub_opt{kidn}(p)).tau0 = tau0;
            NOISE(IndexP_sub_opt{kidn}(p)).dthetadN = dthetadN;
            NOISE(IndexP_sub_opt{kidn}(p)).dRdN = dRdN;
            NOISE(IndexP_sub_opt{kidn}(p)).Fdesign = Fdesign;
            
            % get NEP values smoothed above 10 Hz and below 1 kHz
            % NOISE(kidn).NEPRmin(nT,1) = F,
            % NOISE(kidn).NEPRmin(nT,2) = NEP value
            FR = find(NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,1) < 1e3 & NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,1) > 15);
            [NOISE(IndexP_sub_opt{kidn}(p)).NEPRmin(1,2),indi]      = min(smooth(NEPR(FR),5));
            NOISE(IndexP_sub_opt{kidn}(p)).NEPRmin(1,1)            = NOISE(IndexP_sub_opt{kidn}(p)).NEP(FR(indi),1);
            [NOISE(IndexP_sub_opt{kidn}(p)).NEPthetamin(1,2),indi]  = min(smooth(NEPtheta(FR),5));
            NOISE(IndexP_sub_opt{kidn}(p)).NEPthetamin(1,1)        = NOISE(IndexP_sub_opt{kidn}(p)).NEP(FR(indi),1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            subplot(2,4,1) %Noise
            semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,4),...
                '-','color',kleur(ki,:),'LineWidth',1);%R
            hold on;grid on;xlim([1,1e5]);
            semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,3),...
                '-','color',kleur(ki,:),'LineWidth',2);%Theta
            
            subplot(2,4,2) %Cross PSD
            semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,6),'-',...
                'linewidth',2,'color',kleur(ki,:));hold on;grid on
            semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,5),'-k',...
                'linewidth',1);%fit renormalized
            
            subplot(2,4,3) %Cross PSD
            semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,7),'-',...
                'linewidth',2,'color',kleur(ki,:));hold on;grid on
            
            subplot(2,4,4) %Cross PSD
            semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{1}(:,8),'-',...
                'linewidth',2,'color',kleur(ki,:));hold on;grid on
            
            subplot(2,4,5)% NEP
            loglog(NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,1),         NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,2),'-','Linewidth',1,'color',kleur(ki,:));%R
            grid on;hold on;
            loglog(NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,1),         NOISE(IndexP_sub_opt{kidn}(p)).NEP(:,3),'-','Linewidth',2,'color',kleur(ki,:));%theta
            loglog(NOISE(IndexP_sub_opt{kidn}(p)).NEPRmin(1),       NOISE(IndexP_sub_opt{kidn}(p)).NEPRmin(2),'o','color',kleur(ki,:),'MarkerSize',8);
            loglog(NOISE(IndexP_sub_opt{kidn}(p)).NEPthetamin(1),   NOISE(IndexP_sub_opt{kidn}(p)).NEPthetamin(2),'o','color',kleur(ki,:),'markerfacecolor',kleur(ki,:),'MarkerSize',8);
            
            
            subplot(2,4,6)%lifetime vs P
            plot(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,CrossPSDFit(IndexP_sub_opt{kidn}(p)).tau*1e3,'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            plot(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumin*1e3,'o','color',kleur(ki,:));
            plot(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,CrossPSDFit(IndexP_sub_opt{kidn}(p)).taumax*1e3,'o','color',kleur(ki,:));
            
            subplot(2,4,7)%responsivity            
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,NOISE(IndexP_sub_opt{kidn}(p)).dRdN,'o',...
                'color',kleur(ki,:));hold on;grid on;
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,NOISE(IndexP_sub_opt{kidn}(p)).dthetadN,'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            
            subplot(2,4,8)%NEP vs P
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,NOISE(IndexP_sub_opt{kidn}(p)).NEPRmin(2),'o',...
                'color',kleur(ki,:));hold on;grid on;
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).InternalPower,NOISE(IndexP_sub_opt{kidn}(p)).NEPthetamin(2),'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            MakeGoodFigure(12,9,8)
           
            %Final figure: NEP vs KI and vs Alu
            figure(1234)
            subplot(2,2,1)
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber,NOISE(IndexP_sub_opt{kidn}(p)).NEPthetamin(2),'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            subplot(2,2,2)
            semilogy(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber,NOISE(IndexP_sub_opt{kidn}(p)).NEPRmin(2),'o',...
                'color',kleur(ki,:));hold on;grid on;
            subplot(2,2,3)
            semilogy(AluLength,NOISE(IndexP_sub_opt{kidn}(p)).NEPthetamin(2),'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            subplot(2,2,4)
            semilogy(AluLength,NOISE(IndexP_sub_opt{kidn}(p)).NEPRmin(2),'o',...
                'color',kleur(ki,:));hold on;grid on;
            MakeGoodFigure(9,8,8);
            legc = legc+2;%legend counter
            ki = ki +1; %kleur index
        end %conditional oop if a tau fit is there
        
    end %power loop (p)
    % finding NEP optimum
    % NEPmin.R(kidn) = min([NOISE(IndexP_sub_opt{kidn}).NEPRmin]);%dirty trick:the [] is a single dim array with alternating F(Hz), NEP(W/rt(Hz) etc. Min works as F >> NEP%
    %here we do it nicer
    for mmm = 1: length(IndexP_sub_opt{kidn})
        if ~isempty(NOISE(IndexP_sub_opt{kidn}(mmm)).NEPthetamin)
            Nepjetheta(mmm) = NOISE(IndexP_sub_opt{kidn}(mmm)).NEPthetamin(2);
            NepjeR(mmm) = NOISE(IndexP_sub_opt{kidn}(mmm)).NEPRmin(2);
            tautje(mmm) = CrossPSDFit(IndexP_sub_opt{kidn}(mmm)).tau;
        else
            Nepjetheta(mmm) = NaN;
            NepjeR(mmm) = NaN;
            tautje(mmm) = NaN;
        end
    end
    [NEPmin.theta(kidn),NepIndex] = min(Nepjetheta);
    NEPmin.R(kidn) = min(NepjeR);
    NEPmin.KIDID(kidn) = NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber;
    NEPmin.AluLength(kidn) = AluLength;
    NEPmin.Volume(kidn) = Volume;
    NEPmin.tau(kidn) = tautje(NepIndex);
    clear Nepjetheta NepjeR tautje NepIndex
    
    if notauforthiskid == 0

        figure(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber)
        subplot(2,4,1)
        xlabel('F (Hz)');ylabel('S_x (dBc/Hz)')
        title(['KID_{noise} ' num2str(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber) ', KID_{S21} ' num2str(KID(KIDindresp).KIDnumber)]);
        
        subplot(2,4,2)
        xlabel('F [Hz]');ylabel('S_{Cross}/(dA/dN*d\theta/dN) [1/Hz]');xlim([0.5,1e5]);grid on;
        xlim([10,1e5]);
        title(['Power Dependence @ ' num2str(round(1000*(NOISE(IndexP_sub_opt{kidn}(1)).Temperature))) ' mK'] );
        
        subplot(2,4,3)
        xlabel('F [Hz]');ylabel('S_\theta/(d\theta/dN)^2 [1/Hz]');xlim([0.5,1e5]);grid on;
        xlim([10,1e5]);
        title(['Alu length: ' num2str(AluLength) 'um'])
        
        subplot(2,4,4)
        xlabel('F [Hz]');ylabel('S_A/(dA/dN)^2 [1/Hz]');xlim([0.5,1e5]);grid on;
        xlim([10,1e5]);
        
        subplot(2,4,5)
        legend('NEP_A','NEP_{\theta}','Location','Best');
        xlim([10,1e3]);ylim([0.5e-20 5e-18]);
        xlabel('F (Hz)');ylabel('NEP (W/\surd Hz)')
        title(['V = ' num2str(Volume) '\mum^3, A = ' num2str(Area) '\mum^2']);
        
        subplot(2,4,6)
        xlabel('P  (dBm)');ylabel('\tau (msec)');grid on;ylim([0.01 3]);
        title(['F_{res} = ' num2str(NOISE(IndexP_sub_opt{kidn}(1)).Fres) ' GHz'])
        
        subplot(2,4,7)
        legend('dR/dN','d\theta/dN','Location','Best');
        xlabel('Pi  (dBm)');ylabel('responsivity');grid on;
        
        subplot(2,4,8)%NEP vs P
        legend('NEP_{R}','NEP_{\theta}','Location','NorthWest');
        xlabel('P  (dBm)');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
        
    end
    
    Figfile = [NEPfolder,filesep,...
        'KID',num2str(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber,'%.0f'),'_NEP_Pdep1'];
    MakeGoodFigure(12,9,8,Figfile); %save png and fig
    
    disp(num2str(kidn));
    %for output csv, all params at Max power (pOpt according to noise)
    PoutData(kidn,1) = NOISE(IndexPopt(kidn)).Temperature;
    PoutData(kidn,2) = NOISE(IndexPopt(kidn)).KIDnumber;
    PoutData(kidn,3) = NOISE(IndexPopt(kidn)).InternalPower;
    PoutData(kidn,4) = NOISE(IndexPopt(kidn)).Fres;
    PoutData(kidn,5) = 0.1*round(NOISE(IndexPopt(kidn)).S21min*10);
    PoutData(kidn,6) = round(NOISE(IndexPopt(kidn)).Ql/1e3);
    PoutData(kidn,7) = round(NOISE(IndexPopt(kidn)).Qc/1e3);
    PoutData(kidn,8) = round(NOISE(IndexPopt(kidn)).Qi/1e3);
    PoutData(kidn,9) = Fdesign;
    PoutData(kidn,10) = AluLength;
    PoutData(kidn,11) = Area;
    PoutData(kidn,12) = AluV;
    PoutData(kidn,13) = NOISE(IndexPopt(kidn)).MeanFreqNoise1000Hz;
    PoutData(kidn,14) = dxdN;
    if notauforthiskid == 0
        %NB: the values below are only avilable for some powers, not
        %stored in an index. The raay is initialized with 0's. So
        %min/max works
        PoutData(kidn,15) = min([NOISE(IndexP_sub_opt{kidn}).NEPthetamin]);
        PoutData(kidn,16) = 1e3*max([CrossPSDFit(IndexP_sub_opt{kidn}).tau]);
    end
    %==================================================================================================================================
    %==================================================================================================================================

        
end %KID loop

%Final figure: NEP vs KI and vs Alu
figure(1234)
subplot(2,2,1)
semilogy(NEPmin.KIDID,NEPmin.theta,'xk','MarkerSize',8);  hold on  
xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
title('Phase readout')
subplot(2,2,2)
semilogy(NEPmin.KIDID,NEPmin.R,'+k','MarkerSize',8);  hold on  
xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
title('R readout')
subplot(2,2,3)
semilogy(NEPmin.AluLength,NEPmin.theta,'xk','MarkerSize',8);  hold on 
xlabel('Alu Length  (\mum)');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
title('Phase readout')
subplot(2,2,4)
semilogy(NEPmin.AluLength,NEPmin.R,'+k','MarkerSize',8);  hold on  
xlabel('Alu Length  (\mum)');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
title('R readout')

Figfile = [NEPfolder,filesep,...
    'AllKIDs'];
MakeGoodFigure(9,8,8,Figfile); %save png and fig

save([NEPfolder,filesep,'NEP_P.mat'],'NOISE','IndexP_sub_opt','KIDnumbers','IndexPopt','IndexPsort','CrossPSDNOISE','CrossPSDFit','NEPmin');
PoutHeader = {'T (K)','KIDID','Pint_opt (dBm)','Fres (GHz)','S21 (dB)','Q','Qc','Qi','Fdesign',...
    'AluLength (um)','Area (um^2)','AluV (um^3)','1kHzFnoise','ddxdNqp','NEPtheta','tau (msec)'}';
WriteSRONcsv([NEPfolder,filesep,'NEP_P.csv'],PoutData,PoutHeader,'%.6g')
rmpath([pwd,filesep,'..',filesep,'subroutines']);
