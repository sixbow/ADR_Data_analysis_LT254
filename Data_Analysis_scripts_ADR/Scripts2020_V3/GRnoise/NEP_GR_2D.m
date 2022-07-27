
%

clear variables;close all;clc

%================================================================================
% Input
%================================================================================
ChipInfo_path = ['..' filesep '..' ];; %root path where data is, one higher than the scripts
NEPfolder    =  [filesep 'NEP'];
OptNEPdir    = '\\MARS\kid\KIDonSun\experiments\Entropy ADR\LT197-chip9\2D_BB\2D_BB' ; 
eta_pb = 0.5;
%================================================================================

addpath([pwd,filesep,'..',filesep,'subroutines']);

%load relevant data noise files (level, lifetime etc)
load([ChipInfo_path NEPfolder,filesep,'NEP_2D.mat'],...
    'NOISE','KIDnumbers','IndexP_sub_opt','IndexPopt','CrossPSDFit','NOISEPdeplowT','CrossPSDNOISE','KID');
disp([ 'Low T limit implemented in Cross PSD: ' num2str(NOISEPdeplowT(1).LowTLIM(1)) ' K']);
ChipInfo.path = ChipInfo_path;clear ChipInfo_path; %needed to keep struct ok wrt path names when data is moved in between analysis ruyns

%define NEP folder
NEPfolder = [ChipInfo.path, filesep, 'NEP'] ;
if ~isfolder(NEPfolder)
    mkdir(NEPfolder);
end
%define output array
PoutData = zeros(length(KIDnumbers),17);

for kidn=1:length(KIDnumbers) % LOOP OVER ALL UNIQUE KIDS,
    %construct filename
    legc = 1;%legend counter
    clear lstr
    ki = 1;%index for color plots
    notauforthiskid = 1;
    
   
    for p=1:length(IndexP_sub_opt{kidn})% over Power
        %getting minimum temperature
        %all rest done at Tmin
        [Tmin, Tminind] = min(NOISE(IndexP_sub_opt{kidn}(p)).Temperature); %DO NOT put this outside the P loop, it will cause inconsistencies with the NoiseFinal_GR_2D.m script%
        
        %Getting indexing response file
        KIDID = NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber;
        bla = [KID(:).KIDnumber];
        KIDindresp = (bla == KIDID);
        if sum(KIDindresp) ~= 1
            error('S21 files have multiple readout powers')
        end
        
        Delta = KID(KIDindresp).Delta;%not T dependent in the data!
        Tc = KID(KIDindresp).Tc;%not T dependent in the data!
        Area = KID(KIDindresp).Area;%not T dependent in the data!
        Volume = KID(KIDindresp).Volume;%not T dependent in the data!
        Fres_resp = KID(KIDindresp).Fres(1);%Fres, not T dependent in the data!
        Fres = NOISE(IndexP_sub_opt{kidn}(p)).Fres(1);
        AluLength = KID(KIDindresp).AluLength;
        dxdN = KID(KIDindresp).ddxdNqp(1);

        figure(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber)
        kleur = colormapJetJB(sum( ~isnan(NOISEPdeplowT(kidn).crosslevel) ));
        
        %check if we have a lifetime, otherwise we do not bother
        if ~isnan(NOISEPdeplowT(kidn).crosslevel(p)) %not 0, so we have tau fit and can continue
            notauforthiskid = 0; %tau is there!
            KIDtau_T= 1e-9*KID(KIDindresp).Ql./(pi*KID(KIDindresp).Fres); %tau res vs T in sec
            F       = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,1); %frequency cross PSD
            CPSD    = 10*log10(abs(-1*real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,2)))); %cross PSD
            SR      = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,4);
            Stheta  = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,3);
            %cross fit results
            crosslevel  = CrossPSDFit(IndexP_sub_opt{kidn}(p)).crosslevel(Tminind);
            crosssetup  = CrossPSDFit(IndexP_sub_opt{kidn}(p)).setupnoise(Tminind);
            KIDIDcheck  = NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber;
            if KIDID ~= KIDIDcheck
                error('KID ID problems!!')
            end
            Pint = NOISE(IndexP_sub_opt{kidn}(p)).InternalPower(Tminind);
            Temperature = NOISE(IndexP_sub_opt{kidn}(p)).Temperature(Tminind);
            Fres = NOISE(IndexP_sub_opt{kidn}(p)).Fres(Tminind);
            % kid stuff @ low T
            dthetadN    = NOISEPdeplowT(kidn).dthetadN_Tmin(p)   ;
            dRdN        = NOISEPdeplowT(kidn).dRdN_Tmin(p)        ;
            tau         = NOISEPdeplowT(kidn).Tau(p); %this is the avergae tau at low T
            % put in lowT datastruct som more stuff
            NOISEPdeplowT(kidn).Pint(p)     = Pint;
            NOISEPdeplowT(kidn).tau(p)      = tau;
            NOISEPdeplowT(kidn).T(p)        = Tmin;
            NOISEPdeplowT(kidn).ReadPower(p)= NOISE(IndexP_sub_opt{kidn}(p)).ReadPower(Tminind);
        
            %get reponse and other params @ correct T by interpolation
            if min(KID(KIDindresp).Temperature) <= Tmin && max(KID(KIDindresp).Temperature) > Tmin
                taures =    interp1(KID(KIDindresp).Temperature, KIDtau_T, Tmin);
            elseif min(KID(KIDindresp).Temperature) > Tmin && max(KID(KIDindresp).Temperature) > Tmin
                taures =    KIDtau_T(1);
            else
                error('Repsonse T range terrible, cannot consruct NEP')
            end
           
            % Get the NEP @ min T. NOTE: I do not use the tau_0, Nqp etc
            % from getDARKNEP. It depends on cross PSD and a lifetime fit
            % is probably more reliable
            [NEPR,NEPtheta,NEPcross,Nqp,nqp,NEPGR,tau0,dthetaPdark] = getDARKNEP(dthetadN, dRdN, Stheta, SR, CPSD, F, eta_pb, tau, taures, crosslevel, Volume, Delta, Tc);
            
            % get NEP values smoothed in F range above 10 Hz and below 1 kHz
            FR = find(F < 1e3 & F > 15);
            
            %smooth and find NEP minimum values (p dependent arrays)
            [NOISEPdeplowT(kidn).NEPR(p),indi]      = min(smooth(NEPR(FR),5)); %minimum R NEP
            NOISEPdeplowT(kidn).FNEPR(p)            = F(FR(indi));
            
            [NOISEPdeplowT(kidn).NEPtheta(p),indi]  = min(smooth(NEPtheta(FR),5));
            NOISEPdeplowT(kidn).FNEPtheta(p)        = F(FR(indi));
            NOISEPdeplowT(kidn).Stheta_NEP(p)       = 10.^(Stheta(FR(indi))/20);
            NOISEPdeplowT(kidn).dthetaPdark(p)      = dthetaPdark;
            
            NOISEPdeplowT(kidn).KIDnumber(p)        = NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber; % as check, == kidn

            
            %store NEP curves
            NOISElowT(IndexP_sub_opt{kidn}(p)).NEPF     = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,1); % == F
            NOISElowT(IndexP_sub_opt{kidn}(p)).NEPR     = NEPR;
            NOISElowT(IndexP_sub_opt{kidn}(p)).NEPtheta = NEPtheta;
            NOISElowT(IndexP_sub_opt{kidn}(p)).NEPcross = NEPcross;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lstr{p} = ['P_{read} = ' num2str(NOISEPdeplowT(kidn).ReadPower(p)) ' dBm'];
            subplot(2,3,1) %Noise
            semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,4),...
                '-','color',kleur(ki,:),'LineWidth',1);%R
            hold on;grid on;xlim([1,1e5]);
            semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,3),...
                '-','color',kleur(ki,:),'LineWidth',2);%Theta
            
            subplot(2,3,2) %Cross PSD
            semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,6),'-',...
                'linewidth',2,'color',kleur(ki,:));hold on;grid on
            semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,1),CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{Tminind}(:,5),'-k',...
                'linewidth',1);%fit renormalized
                        
            subplot(2,3,3)%lifetime vs P
            plot(NOISEPdeplowT(kidn).Pint(p),NOISEPdeplowT(kidn).tau(p)*1e3,'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            
            subplot(2,3,4)%NEP vs P
            semilogy(NOISEPdeplowT(kidn).Pint(p),NOISEPdeplowT(kidn).NEPR(p),'o',...
                'color',kleur(ki,:));hold on;grid on;
            semilogy(NOISEPdeplowT(kidn).Pint(p),NOISEPdeplowT(kidn).NEPtheta(p),'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            
            subplot(2,3,5)%resp vs P
            semilogy(NOISEPdeplowT(kidn).Pint(p),NOISEPdeplowT(kidn).dRdN_Tmin(p),'o',...
                'color',kleur(ki,:));hold on;grid on;
            semilogy(NOISEPdeplowT(kidn).Pint(p),NOISEPdeplowT(kidn).dthetadN_Tmin(p),'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            
            subplot(2,3,6)% NEP
            loglog(NOISElowT(IndexP_sub_opt{kidn}(p)).NEPF,	NOISElowT(IndexP_sub_opt{kidn}(p)).NEPR,'-','Linewidth',1,'color',kleur(ki,:));%R
            grid on;hold on;
            loglog(NOISElowT(IndexP_sub_opt{kidn}(p)).NEPF,	NOISElowT(IndexP_sub_opt{kidn}(p)).NEPtheta,'-','Linewidth',2,'color',kleur(ki,:));%theta
            loglog(NOISEPdeplowT(kidn).FNEPR(p),               NOISEPdeplowT(kidn).NEPR(p),'o','color',kleur(ki,:),'MarkerSize',8);
            loglog(NOISEPdeplowT(kidn).FNEPtheta(p),           NOISEPdeplowT(kidn).NEPtheta(p),'o','color',kleur(ki,:),'markerfacecolor',kleur(ki,:),'MarkerSize',8);

            MakeGoodFigure(12,9,8)
           
            %Final figure: NEP vs KI and vs Alu
            figure(1234)
            subplot(2,3,1)%NEP
            semilogy(NOISEPdeplowT(kidn).KIDnumber(p),NOISEPdeplowT(kidn).NEPtheta(p),'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            semilogy(NOISEPdeplowT(kidn).KIDnumber(p),NOISEPdeplowT(kidn).NEPR(p),'s',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            
            subplot(2,3,3)%tau
            semilogy(NOISEPdeplowT(kidn).KIDnumber(p),NOISEPdeplowT(kidn).Tau(p)*1e3,'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            %vs KID L
            subplot(2,3,4)
            semilogy(KID(KIDindresp).AluLength,NOISEPdeplowT(kidn).NEPtheta(p),'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            semilogy(KID(KIDindresp).AluLength,NOISEPdeplowT(kidn).NEPR(p),'s',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            
            subplot(2,3,6)%tau
            semilogy(KID(KIDindresp).AluLength,NOISEPdeplowT(kidn).Tau(p)*1e3,'o',...
                'color',kleur(ki,:),'MarkerFaceColor',kleur(ki,:));hold on;grid on;
            
           MakeGoodFigure(12,9,8);
            
            legc = legc+2;%legend counter
            ki = ki +1; %kleur index
        else %no tau data so no NEP etc
            NOISEPdeplowT(kidn).NEPR(p)     = 1; %minimum R NEP set to very high
            NOISEPdeplowT(kidn).FNEPR(p)	= 0;
            NOISEPdeplowT(kidn).NEPtheta(p) = 1;
            NOISEPdeplowT(kidn).FNEPtheta(p) =0;
            NOISEPdeplowT(kidn).KIDnumber(p) = NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber; % as check, == kidn
            
        end %conditional oop if a tau fit is there
        
    end %power loop (p)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % finish figure    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if notauforthiskid == 0 %do only if data is there

        figure(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber)
        subplot(2,3,1)
        xlabel('F (Hz)');ylabel('S_x (dBc/Hz)');xlim([2,1e5])
        title(['KID_{noise} ' num2str(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber) ', KID_{S21} ' num2str(KID(KIDindresp).KIDnumber) 'KID_{NOISEPdeplowT}' num2str(NOISEPdeplowT(kidn).KIDnumber(p))]);
        
        subplot(2,3,2)
        xlabel('F [Hz]');ylabel('S_{Cross}/(dA/dN*d\theta/dN) [1/Hz]');xlim([0.5,1e5]);grid on;
        xlim([10,1e5]);
        title(['Power Dependence for @T ' num2str(Tmin) ' K'] );
        
        subplot(2,3,3)
        xlabel('Pi  (dBm)');ylabel('\tau (msec)');grid on;ylim([0.01 3]);
        title(['F_{res} = ' num2str(NOISE(IndexP_sub_opt{kidn}(1)).Fres(1)) ' GHz'])

        subplot(2,3,4)%NEP vs P
        legend('NEP_{R}','NEP_{\theta}','Location','NorthWest');
        xlabel('Pi  (dBm)');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
        
        subplot(2,3,5)%NEP vs res
        legend('dR/dN','d\theta/dN','Location','Best');
        xlabel('Pi  (dBm)');ylabel('responsivity');grid on;
        
        subplot(2,3,6)
        legend('NEP_A','NEP_{\theta}','Location','SouthEast');
        xlim([2,1e3]);ylim([0.5e-20 5e-18]);
        xlabel('F (Hz)');ylabel('NEP (W/\surd Hz)')
        title(['V = ' num2str(Volume) '\mum^3, A = ' num2str(Area) '\mum^2, L = ' num2str(AluLength) ' \mum']);

        MakeGoodFigure(12,9,8)
    end
    
    Figfile = [NEPfolder,filesep,...
        'KID',num2str(NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber,'%.0f'),'_NEP_Pdep1'];
    MakeGoodFigure(12,9,8,Figfile); %save png and fig

    disp(num2str(kidn));
    
    if notauforthiskid == 0 %do only if there is any data for at least 1 power
        % Find optimum parameters
        % getting NEP optimum and store in NEPmin
        [NEPmin.theta(kidn),NepIndex]   = min(NOISEPdeplowT(kidn).NEPtheta);
        NEPmin.R(kidn)                  = min(NOISEPdeplowT(kidn).NEPR);
        NEPmin.ReadPower(kidn)          = NOISEPdeplowT(kidn).ReadPower(NepIndex);
        NEPmin.Stheta_NEP(kidn)         = NOISEPdeplowT(kidn).Stheta_NEP(NepIndex);%dthetaPdark
        NEPmin.dthetaPdark(kidn)        = NOISEPdeplowT(kidn).dthetaPdark(NepIndex);%
        NEPmin.KIDID(kidn)              = NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber;%== kidn
        NEPmin.AluLength(kidn)          = AluLength;
        NEPmin.Volume(kidn)             = Volume;
        NEPmin.tau(kidn)                = NOISEPdeplowT(kidn).tau(NepIndex);
        NEPmin.Tmin(kidn)               = Tmin;
        NEPmin.dthetadN(kidn)        	= NOISEPdeplowT(kidn).dthetadN_Tmin(NepIndex);
        NEPmin.dRdN(kidn)               = NOISEPdeplowT(kidn).dRdN_Tmin(NepIndex);
        NEPmin.NEPGRquick(kidn)         = NEPquick2(NEPmin.tau(kidn),Volume);
        NEPmin.dxdN(kidn)               = -1*dxdN; %FIT parameters of the freq. responsivity. 1 is slope dx(Nqp);
        NEPmin.dxdPdark(kidn)           = -1*dxdN * eta_pb * NEPmin.tau(kidn)/Delta; %FIT parameters of the freq. responsivity. 1 is slope dx(Nqp);
        NEPmin.eta_pb                   = eta_pb;
        NEPmin.Delta                    = Delta;


        %for output csv, all params at optimum NEP NepIndex
        PoutData(kidn,1) = Tmin;
        PoutData(kidn,2) = NEPmin.KIDID(kidn);
        PoutData(kidn,3) = NEPmin.ReadPower(kidn);
        PoutData(kidn,4) = NOISE(IndexP_sub_opt{kidn}(NepIndex)).Fres(Tminind);
        PoutData(kidn,5) = round(NOISE(IndexP_sub_opt{kidn}(NepIndex)).S21min(Tminind)*10)/10;
        PoutData(kidn,6) = round(NOISE(IndexP_sub_opt{kidn}(NepIndex)).Ql(Tminind)/1e3);
        PoutData(kidn,7) = round(NOISE(IndexP_sub_opt{kidn}(NepIndex)).Qc(Tminind)/1e3);
        PoutData(kidn,8) = round(NOISE(IndexP_sub_opt{kidn}(NepIndex)).Qi(Tminind)/1e3);
        PoutData(kidn,9) = KID(KIDindresp).Fdesign;
        PoutData(kidn,10) = AluLength;
        PoutData(kidn,11) = Area;
        PoutData(kidn,12) = Volume;
        PoutData(kidn,13) = NOISE(IndexP_sub_opt{kidn}(NepIndex)).MeanFreqNoise1000Hz(Tminind);
        PoutData(kidn,14) = dxdN;
        PoutData(kidn,15) = NEPmin.R(kidn);
        PoutData(kidn,16) = NEPmin.theta(kidn);
        PoutData(kidn,17) = 1e3*NEPmin.tau(kidn);
    else
        NEPmin.theta(kidn)      = 1;
        NEPmin.R(kidn)          = 1;
        NEPmin.ReadPower(kidn)          = NOISEPdeplowT(kidn).ReadPower(NepIndex);
        NEPmin.KIDID(kidn)              = NOISE(IndexP_sub_opt{kidn}(1)).KIDnumber;%== kidn
        NEPmin.AluLength(kidn)          = AluLength;
        NEPmin.Volume(kidn)             = Volume;
        NEPmin.tau(kidn)                = 0;
        NEPmin.Tmin(kidn)               = Tmin;
        NEPmin.dthetadN(kidn)        	= 0;
        NEPmin.dRdN(kidn)               = 0;
        NEPmin.NEPGRquick(kidn)         = 0;
        NEPmin.dxdN(kidn)               = 0; %FIT parameters of the freq. responsivity. 1 is slope dx(Nqp);
        NEPmin.dxdPdark(kidn)           = 0; %FIT parameters of the freq. responsivity. 1 is slope dx(Nqp);

    end
    %==================================================================================================================================
    %==================================================================================================================================
    
    
        
end %KID loop

%Final figure: NEP vs KI and vs Alu
figure(1234)
subplot(2,3,1)
semilogy(NEPmin.KIDID,NEPmin.theta,'ok','MarkerSize',9);  hold on  
semilogy(NEPmin.KIDID,NEPmin.R,'sk','MarkerSize',9); 
xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
legend('Phase readout','R readout')

subplot(2,3,2)%resp
semilogy(NEPmin.KIDID,NEPmin.dthetadN,'or','MarkerSize',8,'MarkerFaceColor','r');  hold on  
semilogy(NEPmin.KIDID,NEPmin.dRdN,'sb','MarkerSize',8,'MarkerFaceColor','b'); 
semilogy([KID.KIDnumber],[KID.dthetadNqpmaxM1],'ok','MarkerSize',8);  hold on  
semilogy([KID.KIDnumber],[KID.dRdNqpmaxM1],'sk','MarkerSize',8); 
xlabel('KID ID');ylabel('responsivity');grid on;
legend('d\theta/dN','dR/dN','max(d\theta/dN) from S21(T)','max(dR/dN) from S21(T)','location','best')

subplot(2,3,3)
semilogy(NEPmin.KIDID,NEPmin.tau*1e3,'ok','MarkerSize',9);  hold on 
xlabel('Alu Length  (\mum)');ylabel('\tau_{qp} (msec.)');grid on;

subplot(2,3,4)
semilogy(NEPmin.AluLength,NEPmin.theta,'ok','MarkerSize',9);  hold on  
semilogy(NEPmin.AluLength,NEPmin.R,'sk','MarkerSize',9); 
xlabel('Alu length (\mum)');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
legend('Phase readout','R readout')

subplot(2,3,5)%AluLength
semilogy(NEPmin.AluLength,NEPmin.dthetadN,'or','MarkerSize',8,'MarkerFaceColor','r');  hold on  
semilogy(NEPmin.AluLength,NEPmin.dRdN,'sb','MarkerSize',8,'MarkerFaceColor','b'); 
semilogy([KID.AluLength],[KID.dthetadNqpmaxM1],'ok','MarkerSize',8);  hold on  
semilogy([KID.AluLength],[KID.dRdNqpmaxM1],'ok','MarkerSize',8); 
xlabel('Alu length (\mum)');ylabel('responsivity');grid on;
legend('d\theta/dN','dR/dN','max(d\theta/dN) from S21(T)','max(dR/dN) from S21(T)','location','best')

subplot(2,3,6)
semilogy(NEPmin.AluLength,NEPmin.tau*1e3,'ok','MarkerSize',9);  hold on 
xlabel('Alu Length  (\mum)');ylabel('\tau_{qp} (msec.)');grid on;ylim([0.1 10])

Figfile = [NEPfolder,filesep,'AllKIDs_DarkNEP'];
MakeGoodFigure(13,9,8,Figfile); %save png and fig

Compare_opticalNEP_DarkNEP(NEPmin,OptNEPdir);
figure(1000)
Figfile = [NEPfolder,filesep,'AllKIDs_DarkOpticalGRNEP'];
MakeGoodFigure(12,9,8,Figfile); %save png and fig

figure(2000)
Figfile = [NEPfolder,filesep,'AllKIDs_DarkOpticalGRNEP_debug'];
MakeGoodFigure(12,9,8,Figfile); %save png and fig

save([NEPfolder,filesep,'NEP_2D.mat'],...
    'NEPmin','NOISElowT','-append');
PoutHeader = {'T (K)','KIDID','Pread (dBm)','Fres (GHz)','S21 (dB)','Q','Qc','Qi','Fdesign',...
    'AluLength (um)','Area (um^2)','AluV (um^3)','1kHzFnoise','ddxdNqp','NEPR','NEPtheta','tau (msec)'}';
WriteSRONcsv([NEPfolder,filesep,'NEP_P.csv'],PoutData,PoutHeader,'%.6g')
rmpath([pwd,filesep,'..',filesep,'subroutines']);
