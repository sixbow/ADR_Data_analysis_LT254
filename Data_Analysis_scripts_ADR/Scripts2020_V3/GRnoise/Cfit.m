classdef Cfit < handle
    %Cfit is the base class for our data analysis.
    %
    
    properties
        %Data vars
        NOISE
        CrossPSDFit
        Cross_tau
        Cross_taumin
        Cross_taumax
        IndexPsort
        IndexPopt
        freq
        % matrices of generated data.
        Sff
        Sff_sys_noise
        Sff_minusTLS
        fRing
        %Looping vars
        kidn_iter
        power_iter
        Tbath_iter
        %Fitting vars
        TLSfit
        CBfit
        fguess
        resnormthreshold
        Goodfit
        Fknee
        Sknee
        % Anonymous functions
        fTLS
        fCB
        fCBmintau
        fCBmaxtau
        fline
        fTotal
        fTLS_shallow
        fTLS_filled
        fCB_filled
        %Plotting vars
        fig % Figure handle.
        ax  % Axes handle.]
        test 
    end
    
    methods
        function obj = Cfit(FFTsubsubdir,filename,filenameCPSDfit,TLSfitmin,TLSfitmax,CBfitmin,CMfitmax)
            %Cfit Construct an instance of this class
            % First argument gives the path of the data from 2 folders up. 
            % Second argument is the type of data that you want to load.
            % 'NOISE_P.mat' or 'NOISE_2D.mat'
            addpath('..\subroutines')
            ChipInfo_path = ['..' filesep '..' filesep ]; %root path where data is, one higher than the scripts       
            obj.NOISE = load([ChipInfo_path FFTsubsubdir filesep filename],'NOISE').NOISE;
            obj.IndexPsort = load([ChipInfo_path FFTsubsubdir filesep filename],'IndexPsort').IndexPsort;
            obj.IndexPopt = load([ChipInfo_path FFTsubsubdir filesep filename],'IndexPopt').IndexPopt;
            obj.CrossPSDFit = load([ChipInfo_path FFTsubsubdir filesep filenameCPSDfit],'CrossPSDFit').CrossPSDFit;
            obj.freq = obj.NOISE(1).FFTnoise{1,1}(:,1);
            obj.kidn_iter = 1:6; % Default value
            obj.power_iter = obj.IndexPsort{obj.kidn_iter,1};
            obj.Tbath_iter = 1:14;
            %-----TLS fitting parameters (Default!)
            obj.TLSfit.min = TLSfitmin;% [Hz]
            obj.TLSfit.max= TLSfitmax;% [Hz]
            [~,obj.TLSfit.mini] = min(abs(obj.freq-obj.TLSfit.min)); % This gives the min index
            [~,obj.TLSfit.maxi] = min(abs(obj.freq-obj.TLSfit.max)); % This gives the max index
            obj.TLSfit.C0 = [(10^-15) 1.5];% [CTLS gamma]  linear scale. CTLS at 1 Hz , power with what  the power law decreases. gamma TLS ~0.5
 			%-------/End:TLS fitting parameters (Default!)-----------------
            %+++++|Begin:Combined fitting parameters (Default!)------------
            obj.CBfit.min = CBfitmin;% [Hz]
            obj.CBfit.max = CMfitmax;%  [Hz]
            [~,obj.CBfit.mini] = min(abs(obj.freq-obj.CBfit.min)); % "
            [~,obj.CBfit.maxi] = min(abs(obj.freq-obj.CBfit.max)); % "
            obj.CBfit.C0 = [(10^-18) 1.591549430918954e-06];% CTLS at 1 Hz , power with what  the power law decreases. TLS ~0.5
            obj.CBfit.lb = [10^-20 (1/(2*pi*obj.CBfit.max))];
            obj.CBfit.ub = [10^-15 (1/(2*pi*obj.CBfit.min))];
            obj.fguess = [1 100000]; % [Hz] initial guess for the intersection point of the TLS and GR noise. 
			%-------/End:Combined fitting parameters (Default!)------------
            obj.resnormthreshold = 2000;


            disp('Initialized Cfit object!')
        end % End constructor
        %+++++:Begin fitting and calcuations.-------------------------------
        
        
        function  [Ctls,gamma]= genfitsingle(obj,kidn,Pindex,nT)
            %genfit method fits the TLS fit and the GR fits 
            p = obj.findp(kidn,Pindex);
            obj.Sff{kidn,Pindex,nT} = obj.NOISE(p).FFTnoise{nT}(:,4);
            setup_noise_eval_freq = 350000;% [Hz]
            [~,Snoise_i] = min(abs(obj.freq-setup_noise_eval_freq)); % "
            obj.Sff_sys_noise{kidn,Pindex,nT} = obj.NOISE(p).FFTnoise{nT}(Snoise_i,4); % This is the setup noise level evaluated at 350000Hz
            %Updating and Setting anonymous functions
            obj.UpA(kidn,Pindex,nT)
            %disp('hoi')
			%+++++|Begin:TLS fitting---------------------------------------
            dBlogDataX = log10(obj.freq);
            dBlogDataY = 10.*log10(obj.NOISE(p).FFTnoise{nT}(:,4));
            TLS_coof_dBlog{kidn,Pindex,nT} = polyfit(dBlogDataX(obj.TLSfit.mini:obj.TLSfit.maxi),dBlogDataY(obj.TLSfit.mini:obj.TLSfit.maxi),1);
            obj.TLSfit.Ctls{kidn,Pindex,nT} =  power(10,(TLS_coof_dBlog{kidn,Pindex,nT}(2)/10));
            obj.TLSfit.gamma{kidn,Pindex,nT} = -(TLS_coof_dBlog{kidn,Pindex,nT}(1)/10);
            obj.TLSfit.C{kidn,Pindex,nT} = [obj.TLSfit.Ctls{kidn,Pindex,nT} obj.TLSfit.gamma{kidn,Pindex,nT}];
            %disp(string(obj.TLSfit.C{kidn,Pindex,nT}(1)))
            Ctls = obj.TLSfit.Ctls;% Output to user
            gamma = obj.TLSfit.gamma;% Output to user
            
            %disp('doei')
			%-------/End:TLS fitting---------------------------------------
         
			%+++++|Begin:Combined fitting----------------------------------
			% Fit GR-noise + TLS ~400Hz - 100KHz
            %disp(string(obj.TLSfit.C{kidn,Pindex,nT}(1)))
            
            %disp(string(obj.TLSfit.C{kidn,Pindex,nT}(2)))
            %disp(obj.fTLS)
            obj.Sff_minusTLS{kidn,Pindex,nT} = obj.Sff{kidn,Pindex,nT}-obj.fTLS(obj.TLSfit.C{kidn,Pindex,nT},obj.freq); % subtract the TLS line..
            
            %Fitting GR noise spectrum 
            [obj.CBfit.C{kidn,Pindex,nT},obj.CBfit.resnorm{kidn,Pindex,nT}] = LLS_CB_SdB(obj.freq(obj.CBfit.mini:obj.CBfit.maxi),obj.Sff_minusTLS{kidn,Pindex,nT}(obj.CBfit.mini:obj.CBfit.maxi),obj.fCB,obj.CBfit.C0,obj.CBfit.lb,obj.CBfit.ub);
            obj.CBfit.Cgr{kidn,Pindex,nT} = obj.CBfit.C{kidn,Pindex,nT}(1);
            obj.CBfit.Tauqp{kidn,Pindex,nT} = obj.CBfit.C{kidn,Pindex,nT}(2);
            disp(fprintf('genfitsingle - ID(%i,%i,%i): resnorm=%2.4e',kidn,Pindex,nT,obj.CBfit.resnorm{kidn,Pindex,nT}));
            %disp('Fitted GR!')    
            
            
            
            
            
            
            %-------/End:Combined fitting---------------------------------- 
        end% End genfitsingle
        
        function  [Ctls,gamma]= genfitCrosssingle(obj,kidn,Pindex,nT)
            %genfit method fits the TLS fit and the GR fits 
            p = obj.findp(kidn,Pindex);
            obj.Sff{kidn,Pindex,nT} = obj.NOISE(p).FFTnoise{nT}(:,4);
            setup_noise_eval_freq = 350000;% [Hz]
            [~,Snoise_i] = min(abs(obj.freq-setup_noise_eval_freq)); % "
            obj.Sff_sys_noise{kidn,Pindex,nT} = obj.NOISE(p).FFTnoise{nT}(Snoise_i,4); % This is the setup noise level evaluated at 350000Hz
            %Updating and Setting anonymous functions
            obj.UpCross(kidn,Pindex,nT)
            %disp('hoi')
			%+++++|Begin:TLS fitting---------------------------------------
            dBlogDataX = log10(obj.freq);
            dBlogDataY = 10.*log10(obj.NOISE(p).FFTnoise{nT}(:,4));
            TLS_coof_dBlog{kidn,Pindex,nT} = polyfit(dBlogDataX(obj.TLSfit.mini:obj.TLSfit.maxi),dBlogDataY(obj.TLSfit.mini:obj.TLSfit.maxi),1);
            obj.TLSfit.Ctls{kidn,Pindex,nT} =  power(10,(TLS_coof_dBlog{kidn,Pindex,nT}(2)/10));
            obj.TLSfit.gamma{kidn,Pindex,nT} = -(TLS_coof_dBlog{kidn,Pindex,nT}(1)/10);
            obj.TLSfit.C{kidn,Pindex,nT} = [obj.TLSfit.Ctls{kidn,Pindex,nT} obj.TLSfit.gamma{kidn,Pindex,nT}];
            %disp(string(obj.TLSfit.C{kidn,Pindex,nT}(1)))
            Ctls = obj.TLSfit.Ctls;% Output to user
            gamma = obj.TLSfit.gamma;% Output to user
            
            %disp('doei')
			%-------/End:TLS fitting---------------------------------------
         
			%+++++|Begin:Combined fitting----------------------------------
			% Fit GR-noise + TLS ~400Hz - 100KHz
            %disp(string(obj.TLSfit.C{kidn,Pindex,nT}(1)))
            
            %disp(string(obj.TLSfit.C{kidn,Pindex,nT}(2)))
            %disp(obj.fTLS)
            obj.Sff_minusTLS{kidn,Pindex,nT} = obj.Sff{kidn,Pindex,nT}-obj.fTLS(obj.TLSfit.C{kidn,Pindex,nT},obj.freq); % subtract the TLS line..
            fit=0;
            fail = 0;
            %Fitting GR noise spectrum 
            try
            [obj.CBfit.C{kidn,Pindex,nT},obj.CBfit.resnorm{kidn,Pindex,nT}] = LLS_CB_SdB_single(obj.freq(obj.CBfit.mini:obj.CBfit.maxi),obj.Sff_minusTLS{kidn,Pindex,nT}(obj.CBfit.mini:obj.CBfit.maxi),obj.fCB,obj.CBfit.C0(1),obj.CBfit.lb(1),obj.CBfit.ub(1));
            obj.CBfit.Cgr{kidn,Pindex,nT} = obj.CBfit.C{kidn,Pindex,nT}(1);
            obj.CBfit.Tauqp{kidn,Pindex,nT} = obj.Cross_tau{kidn,Pindex,nT};
            disp(fprintf('genfitsingle - ID(%i,%i,%i): resnorm=%2.4e',kidn,Pindex,nT,obj.CBfit.resnorm{kidn,Pindex,nT}));
            fit = fit+1;
            catch e
                fprintf(1,'Fit waarschijnlijk onmogelijk want Tau = Nan -> Setting result to NaN also Message was:\n%s',e.message);
                obj.CBfit.C{kidn,Pindex,nT} = NaN;
                obj.CBfit.resnorm{kidn,Pindex,nT} = NaN;
                obj.CBfit.Cgr{kidn,Pindex,nT} = NaN;
                obj.CBfit.Tauqp{kidn,Pindex,nT} = NaN;
            fail = fail +1;
            end
            disp(['fit :' string(fit)])
            disp(['fail :' string(fail)])
            disp(['fitfail ratio :' string(fit/(fail+fit))])
            
            %disp('Fitted GR!')    
            
            
            
            
            
            
            %-------/End:Combined fitting---------------------------------- 
        end% End genfitCrosssingle
        
        function UpA(obj,kidn,Pindex,nT)
            obj.fTLS = @(C_v,fdata)C_v(1).*power(fdata,-1*C_v(2)); %Model we use to fit C_v is the constants vector that we are trying to find.
            obj.fCB  = @(C_v,fdata) C_v(1)./((1+power((2.*pi.*fdata.*C_v(2)),2)))+ obj.Sff_sys_noise{kidn,Pindex,nT};
            obj.fline = @(C_v,xdata)C_v(1).*xdata+C_v(2);
            obj.fTotal = @(C_v_TLS,C_v_CB,fdata)obj.fTLS(C_v_TLS,fdata)+obj.fCB(C_v_CB,fdata);
        end
        function UpCross(obj,kidn,Pindex,Tindex)% This uses the data from the cross script.
            obj.fTLS = @(C_v,fdata)C_v(1).*power(fdata,-1*C_v(2)); %Model we use to fit C_v is the constants vector that we are trying to find.
            obj.fCB  = @(C_v,fdata) C_v(1)./((1+power((2.*pi.*fdata.*obj.Cross_tau{kidn,Pindex,Tindex}),2)))+ obj.Sff_sys_noise{kidn,Pindex,Tindex};
            obj.fCBmintau  = @(C_v,fdata) C_v(1)./((1+power((2.*pi.*fdata.*obj.Cross_taumin{kidn,Pindex,Tindex}),2)))+ obj.Sff_sys_noise{kidn,Pindex,Tindex};
            obj.fCBmaxtau  = @(C_v,fdata) C_v(1)./((1+power((2.*pi.*fdata.*obj.Cross_taumax{kidn,Pindex,Tindex}),2)))+ obj.Sff_sys_noise{kidn,Pindex,Tindex};
            obj.fline = @(C_v,xdata)C_v(1).*xdata+C_v(2);
            obj.fTotal = @(C_v_TLS,C_v_CB,fdata)obj.fTLS(C_v_TLS,fdata)+obj.fCB(C_v_CB,fdata);
            obj.fTLS_shallow = @(fdata) (10^-16).*power(fdata,-0.3);
        end
        function UpA_filled(obj,kidn,Pindex,nT)
            a = obj.TLSfit.C{kidn,Pindex,nT}(1);
            b = obj.TLSfit.C{kidn,Pindex,nT}(2);
            c = obj.CBfit.C{kidn,Pindex,nT}(1);
            d = obj.CBfit.C{kidn,Pindex,nT}(2);
            
            obj.fTLS_filled = @(fdata)a.*power(abs(fdata),-1*b);
            obj.fCB_filled = @(fdata)c./((1+power((2.*pi.*fdata.*d),2)))+ obj.Sff_sys_noise{kidn,Pindex,nT};
        end
        function genCross_tau(obj,kidn_iter,Pindex_iter,T_iter) % With this you can generate the entries as needed
            for kidn=kidn_iter
                for Pindex = Pindex_iter
                    for Tindex = T_iter
                        p = obj.findp(kidn,Pindex);
                        obj.Cross_tau{kidn,Pindex,Tindex} = obj.CrossPSDFit(p).tau(Tindex);
                        obj.Cross_taumin{kidn,Pindex,Tindex} = obj.CrossPSDFit(p).taumin(Tindex);
                        obj.Cross_taumax{kidn,Pindex,Tindex}  = obj.CrossPSDFit(p).taumax(Tindex);
                    end
                end
            end 
        end
        function genfRing(obj,kidn,Pindex,nT)
        p = obj.findp(kidn,Pindex);
        obj.fRing{kidn,Pindex,nT} = ((obj.NOISE(p).Fres(nT)*10^9)/(2*obj.NOISE(p).Ql(nT)));
        end

        function genFknee(obj,kidn,Pindex,nT)
            %genFknee method finds the Fknee based on the fits.
            % This is done such that fits that are bad are rejected. Also
            % Also if Fknee is greater than rolloff of Tau qp it is
            % rejected as Fknee.
            if obj.CBfit.resnorm{kidn,Pindex,nT} < obj.resnormthreshold
            obj.Fknee{kidn,Pindex,nT} = power(obj.TLSfit.Ctls{kidn,Pindex,nT}/obj.CBfit.Cgr{kidn,Pindex,nT},(1/obj.TLSfit.gamma{kidn,Pindex,nT}));
            obj.Sknee{kidn,Pindex,nT} = obj.CBfit.Cgr{kidn,Pindex,nT};
            obj.Goodfit{kidn,Pindex,nT} = 1;
                if obj.Fknee{kidn,Pindex,nT}> (1/(2*pi*obj.CBfit.Tauqp{kidn,Pindex,nT}))
                obj.Fknee{kidn,Pindex,nT} = NaN;
                obj.Sknee{kidn,Pindex,nT} = NaN;
                obj.Goodfit{kidn,Pindex,nT} = 0;
                end
            else
            obj.Fknee{kidn,Pindex,nT} = NaN;
            obj.Sknee{kidn,Pindex,nT} = NaN;
            obj.Goodfit{kidn,Pindex,nT} = 0;
            end
        end
        %-----:End fitting and calcuations.---------------------------------
        %>>>>>:Begin Plotting fuctions.-------------------------------------
        function  obj = init_figax(obj,fig_n,ax_n,option)
            obj.fig(fig_n) = figure;
            obj.ax(ax_n) = gca;
            if strcmp(option,'linlin')
                set(obj.ax(ax_n),'XScale','linear','YScale','linear')
            elseif strcmp(option,'linlog')
                set(obj.ax(ax_n),'XScale','linear','YScale','log')
            elseif strcmp(option,'loglin')
                set(obj.ax(ax_n),'XScale','log','YScale','linear')
            elseif strcmp(option,'loglog')
                set(obj.ax(ax_n),'XScale','log','YScale','log')
            else
                disp('init_figax: I am confused!')
            end
            hold(obj.ax(ax_n),'on') 
        end
        
        
        
        function [figh_out,axh_out] = plotsingle(obj,fig_n,ax_n,kidn,Pindex,nT,SW,marker,Colorcell,handleVisible)
            %plots to figure currently in focus.
            p = obj.findp(kidn,Pindex);
            if SW.CrossTau
                UpCross(obj,kidn,Pindex,nT)
            else
                UpA(obj,kidn,Pindex,nT)
            end
            obj.Sff{kidn,Pindex,nT} = obj.NOISE(p).FFTnoise{nT}(:,4);
            disp(fprintf('plotsingle:Plotting data with p=%i and nT=%i',p,nT));
            disp(fprintf('genfitsingle - ID(%i,%i,%i): resnorm=%2.4e',kidn,Pindex,nT,obj.CBfit.resnorm{kidn,Pindex,nT}));
            toplot = obj.Sff{kidn,Pindex,nT} > 0;
            figure(obj.fig(fig_n))
            if SW.plotdata
            plot(obj.ax(ax_n),obj.freq(toplot),lintodb(obj.Sff{kidn,Pindex,nT}(toplot)),'-o','MarkerFaceColor',Colorcell{1},'LineWidth',3,'MarkerSize',3,'Color',Colorcell{1},'HandleVisibility',handleVisible{1});
            end
            %plot()
            xlabel('f[Hz]','Interpreter','latex','FontSize',30);ylabel('$S_{F}/F^{2} [dBc/Hz]$','Interpreter','latex','FontSize',30)
            %Old: %xlim([0.5,1e5]);grid on;ylim([-220,-140])
            xlim([0.1,5e5]);grid on;ylim([-200,-140])
            %title(append('KID#',string(obj.NOISE(p).KIDnumber)," |Power ",string(obj.NOISE(p).ReadPower),"dBm")); 
            %+++++|Begin:Add TLS line--------------------------------------
            if SW.plottls
            %disp(string(size(obj.TLSfit.C)))
            %disp(string(obj.TLSfit.C{kidn,Pindex,nT}(1)))
            plot(obj.ax(ax_n),obj.freq(toplot),lintodb(obj.fTLS(obj.TLSfit.C{kidn,Pindex,nT},obj.freq(toplot))),'-.','Color','red','LineWidth',3,'HandleVisibility',handleVisible{2});
            end
            %-------/End:Add TLS line--------------------------------------
            %+++++|Begin:Plot CB fit---------------------------------------
			if SW.plotgr
            plot(obj.ax(ax_n),obj.freq(toplot),lintodb(obj.fCB(obj.CBfit.C{kidn,Pindex,nT},obj.freq(toplot))),':','Color',Colorcell{3},'Color','green','LineWidth',3,'HandleVisibility',handleVisible{3})
            end
            if SW.plottauminmax
            plot(obj.ax(ax_n),obj.freq(toplot),lintodb(obj.fCBmintau(obj.CBfit.C{kidn,Pindex,nT},obj.freq(toplot))),':','Color',Colorcell{4},'LineWidth',0.8,'HandleVisibility','off')
            plot(obj.ax(ax_n),obj.freq(toplot),lintodb(obj.fCBmaxtau(obj.CBfit.C{kidn,Pindex,nT},obj.freq(toplot))),':','Color',Colorcell{4},'LineWidth',0.8,'HandleVisibility','off')
            end
            
            disp('De kleur is:')
            disp(Colorcell{3});
            %-------/End:Plot CB fit---------------------------------------
            %+++++|Begin:Total line----------------------------------------
			if SW.plotTotal
            plot(obj.ax(ax_n),obj.freq(toplot),lintodb(obj.fTotal(obj.TLSfit.C{kidn,Pindex,nT},obj.CBfit.C{kidn,Pindex,nT},obj.freq(toplot))),'--','LineWidth',2,'Color',Colorcell{4},'HandleVisibility',handleVisible{4})
            end
            %-------/End:Total line----------------------------------------
            %+++++|Begin:Total+TSL shallow-------------------------------------
            if SW.plotTotalplusTLSshallow
            plot(obj.ax(ax_n),obj.freq(toplot),lintodb((obj.fTotal(obj.TLSfit.C{kidn,Pindex,nT},obj.CBfit.C{kidn,Pindex,nT},obj.freq(toplot))+obj.fTLS_shallow(obj.freq(toplot)))./(1+(obj.freq(toplot)./obj.fRing{kidn,Pindex,nT}).^2)),'--','LineWidth',4,'Color','m','HandleVisibility',handleVisible{4})
            end
            if SW.plotTotalplusTLSshallow
            plot(obj.ax(ax_n),obj.freq(toplot),lintodb(obj.fTLS_shallow(obj.freq(toplot))./(1+(obj.freq(toplot)./obj.fRing{kidn,Pindex,nT}).^2)),'.-','LineWidth',1,'Color','m','HandleVisibility',handleVisible{4})
            end
            %-------/End:Total+TSL shallow-------------------------------------

            
            
            
            
            %+++++|Begin:Visualize Fknee-----------------------------------
			if SW.plotFknee
            plot(obj.ax(ax_n),obj.Fknee{kidn,Pindex,nT},lintodb(obj.Sknee{kidn,Pindex,nT}),marker,'Color',Colorcell{5},'MarkerSize',8,'LineWidth',2,'HandleVisibility',handleVisible{5})
            end
            %-------/End:Visualize Fknee-----------------------------------
            figh_out = obj.fig(fig_n)
            axh_out = obj.ax(ax_n)
            
            %+++++|Begin:Line ring time------------------------------------
			if SW.plotringline
                xline(obj.ax(ax_n),obj.fRing{kidn,Pindex,nT},'--','HandleVisibility',handleVisible{6});
            end

            %-------/End:Line ring time------------------------------------

            %+++++|Begin:fitting ranges------------------------------------
            if SW.cyanfit
			xline(obj.freq(obj.TLSfit.mini),'--','Color','c','LineWidth',1,'HandleVisibility','off')
            xline(obj.freq(obj.TLSfit.maxi),'--','Color','c','LineWidth',1,'HandleVisibility','off')
            end
            if SW.magentafit
            xline(obj.freq(obj.CBfit.mini),'--','Color','m','LineWidth',1,'HandleVisibility','off')
            xline(obj.freq(obj.CBfit.maxi),'--','Color','m','LineWidth',1,'HandleVisibility','off')
            end
       
            
            
            %-------/End:fitting ranges------------------------------------ 
        end
        
        function plottempfknee(obj,fig_n,ax_n,kidn,Pindex)
        figure(obj.fig(fig_n))
        for Tindex = 1:14
        Tvec(Tindex) = obj.getT(kidn,Pindex,Tindex);
        Fknee(Tindex) = obj.Fknee{kidn,Pindex,Tindex};
        end
        plot(obj.ax(ax_n),Tvec,Fknee,'LineWidth',2)
        end
        
        function plotQifknee(obj,fig_n,ax_n,kidn,Pindex)
        figure(obj.fig(fig_n))
        for Tindex = 1:14
        Qivec(Tindex) = obj.getQi(kidn,Pindex,Tindex);
        Fknee(Tindex) = obj.Fknee{kidn,Pindex,Tindex};
        end
        plot(obj.ax(ax_n),Qivec,Fknee,'LineWidth',2)
        end
        
        function plottautemp(obj,fig_n,ax_n,kidn,Pindex,Tcolors,markersize)
        figure(obj.fig(fig_n))
            for Tindex = 1:14 
            Tvec(Tindex) = obj.getT(kidn,Pindex,Tindex);
            Crosstauvec(Tindex) =  (obj.Cross_tau{kidn,Pindex,Tindex}*1000);
            Crosstauminvec(Tindex) =  (obj.Cross_taumin{kidn,Pindex,Tindex}*1000);
            Crosstaumaxvec(Tindex) =  (obj.Cross_taumax{kidn,Pindex,Tindex}*1000);
            tauvec_own(Tindex) = (obj.CBfit.Tauqp{kidn,Pindex,Tindex}*1000);% to convert to milliseconds.
            plot(obj.ax(ax_n),Tvec(Tindex),tauvec_own(Tindex),'x','LineWidth',1,'MarkerFaceColor',Tcolors(Tindex,:),'MarkerSize',markersize,'Color',Tcolors(Tindex,:))
            plot(obj.ax(ax_n),Tvec(Tindex),Crosstauvec(Tindex),'o','LineWidth',1,'MarkerFaceColor',Tcolors(Tindex,:),'MarkerSize',markersize,'Color',Tcolors(Tindex,:))
            plot(obj.ax(ax_n),Tvec(Tindex),Crosstauminvec(Tindex),'o','LineWidth',1,'MarkerSize',markersize,'Color',Tcolors(Tindex,:))
            plot(obj.ax(ax_n),Tvec(Tindex),Crosstaumaxvec(Tindex),'o','LineWidth',1,'MarkerSize',markersize,'Color',Tcolors(Tindex,:))
            end
        end
       
        function plotcomptau(obj)
            %
            
        end
        function plotmulti(obj)
            %
            
        end
        %<<<<<End Plotting fuctions.---------------------------------------
		
        %+++++|Begin:Getters-----------------------------------------------
        function Temp = getT(obj,kidn,Pindex,nT)
            p = obj.findp(kidn,Pindex);
            Temp = obj.NOISE(p).Temperature(nT);
        end
        function Pread = getPread(obj,kidn,Pindex,~)
            p = obj.findp(kidn,Pindex);
            Pread = obj.NOISE(p).ReadPower;
        end
        function Pread = getFread(obj,kidn,Pindex,~)
            p = obj.findp(kidn,Pindex);
            Pread = obj.NOISE(p).Fread;
        end
        function Pint = getPint(obj,kidn,Pindex,nT)
            p = obj.findp(kidn,Pindex);
            Pint = obj.NOISE(p).InternalPower(nT);
        end
        function out = getQi(obj,kidn,Pindex,nT)
            p = obj.findp(kidn,Pindex);
            out = obj.NOISE(p).Qi(nT);
        end
        function out = getQc(obj,kidn,Pindex,nT)
            p = obj.findp(kidn,Pindex);
            out = obj.NOISE(p).Qc(nT);
        end
        function out = getFres(obj,kidn,Pindex,nT)
            p = obj.findp(kidn,Pindex);
            out = obj.NOISE(p).Fres(nT);
        end
        function out = getTauRes(obj,kidn,Pindex,nT)
            p = obj.findp(kidn,Pindex);
            out = obj.NOISE(p).TauRes(nT);
        end
        function out = getS21min(obj,kidn,Pindex,nT)
            p = obj.findp(kidn,Pindex);
            out = obj.NOISE(p).S21min(nT);
        end
		%-------/End:Getters-----------------------------------------------
		
        
        %+++++|Begin:Setters-----------------------------------------------
		%-------/End:Setters-----------------------------------------------
        
        %+++++|Begin:Auxillary func----------------------------------------
        function p = findp(obj,kidn,Pindex)% this tells you what p belongs to what power index.
            p = obj.IndexPsort{kidn,1}(Pindex);
        end
        
        function Pindex = findPiopt(obj,kidn)% this tells you what p belongs to what power index.
            for i = 1:6
            Pindex(i) = find(obj.findp(i,1:7)==obj.IndexPopt(i));
            end
            Pindex = Pindex(kidn);
        end
        
		%-------/End:Auxillary func----------------------------------------

        
        
        
        
    end
end

