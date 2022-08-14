classdef Cfit < handle
    %Cfit is the base class for our data analysis.
    %
    
    properties
        %Data vars
        NOISE
        IndexPsort
        IndexPopt
        freq
        % matrices of generated data.
        Sff
        Sff_sys_noise
        Sff_minusTLS
        %Looping vars
        kidn_iter
        power_iter
        Tbath_iter
        %Fitting vars
        TLSfit
        CBfit
        fguess
        Fknee
        Sknee
        % Anonymous functions
        fTLS
        fCB
        fline
        fTLS_filled
        fCB_filled
        %Plotting vars
        fig % Figure handle.
        ax  % Axes handle.]
        test 
    end
    
    methods
        function obj = Cfit(FFTsubsubdir,filename)
            %Cfit Construct an instance of this class
            % First argument gives the path of the data from 2 folders up. 
            % Second argument is the type of data that you want to load.
            % 'NOISE_P.mat' or 'NOISE_2D.mat'
            addpath('..\subroutines')
            ChipInfo_path = ['..' filesep '..' filesep ]; %root path where data is, one higher than the scripts       
            obj.NOISE = load([ChipInfo_path FFTsubsubdir filesep filename],'NOISE').NOISE;
            obj.IndexPsort = load([ChipInfo_path FFTsubsubdir filesep filename],'IndexPsort').IndexPsort;
            obj.IndexPopt = load([ChipInfo_path FFTsubsubdir filesep filename],'IndexPopt').IndexPopt;
            
            obj.freq = obj.NOISE(1).FFTnoise{1,1}(:,1);
            obj.kidn_iter = 1:6; % Default value
            obj.power_iter = obj.IndexPsort{obj.kidn_iter,1};
            obj.Tbath_iter = 1:14;
            %-----TLS fitting parameters (Default!)
            obj.TLSfit.min = 0.8;% [Hz]
            obj.TLSfit.max= 30;% [Hz]
            [~,obj.TLSfit.mini] = min(abs(obj.freq-obj.TLSfit.min)); % This gives the min index
            [~,obj.TLSfit.maxi] = min(abs(obj.freq-obj.TLSfit.max)); % This gives the max index
            obj.TLSfit.C0 = [(10^-15) 1.5];% [CTLS gamma]  linear scale. CTLS at 1 Hz , power with what  the power law decreases. gamma TLS ~0.5
 			%-------/End:TLS fitting parameters (Default!)-----------------
            %+++++|Begin:Combined fitting parameters (Default!)------------
            obj.CBfit.min = 1000;% [Hz]
            obj.CBfit.max = 150000;%  [Hz]
            [~,obj.CBfit.mini] = min(abs(obj.freq-obj.CBfit.min)); % "
            [~,obj.CBfit.maxi] = min(abs(obj.freq-obj.CBfit.max)); % "
            obj.CBfit.C0 = [(10^-18) 1.591549430918954e-06];% CTLS at 1 Hz , power with what  the power law decreases. TLS ~0.5
            obj.CBfit.lb = [10^-20 (1/(2*pi*obj.CBfit.max))];
            obj.CBfit.ub = [10^-15 (1/(2*pi*obj.CBfit.min))];
            obj.fguess = [1 100000]; % [Hz] initial guess for the intersection point of the TLS and GR noise. 
			%-------/End:Combined fitting parameters (Default!)------------
            obj.test = 3.14;
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
            obj.CBfit.C{kidn,Pindex,nT} = LLS_CB_SdB(obj.freq(obj.CBfit.mini:obj.CBfit.maxi),obj.Sff_minusTLS{kidn,Pindex,nT}(obj.CBfit.mini:obj.CBfit.maxi),obj.fCB,obj.CBfit.C0,obj.CBfit.lb,obj.CBfit.ub);
            obj.CBfit.Cgr{kidn,Pindex,nT} = obj.CBfit.C{kidn,Pindex,nT}(1);
            obj.CBfit.Tauqp{kidn,Pindex,nT} = obj.CBfit.C{kidn,Pindex,nT}(2);
            
            disp('Fitted GR!')    
            
            
            
            
            
            
            %-------/End:Combined fitting---------------------------------- 
        end% End genfitsingle
        
        function UpA(obj,kidn,Pindex,nT)
            obj.fTLS = @(C_v,fdata)C_v(1).*power(fdata,-1*C_v(2)); %Model we use to fit C_v is the constants vector that we are trying to find.
            obj.fCB  = @(C_v,fdata) C_v(1)./((1+power((2.*pi.*fdata.*C_v(2)),2)))+ obj.Sff_sys_noise{kidn,Pindex,nT};
            obj.fline = @(C_v,xdata)C_v(1).*xdata+C_v(2);
        end
        function UpA_filled(obj,kidn,Pindex,nT)
            a = obj.TLSfit.C{kidn,Pindex,nT}(1);
            b = obj.TLSfit.C{kidn,Pindex,nT}(2);
            c = obj.CBfit.C{kidn,Pindex,nT}(1);
            d = obj.CBfit.C{kidn,Pindex,nT}(2);
            
            obj.fTLS_filled = @(fdata)a.*power(abs(fdata),-1*b);
            obj.fCB_filled = @(fdata)c./((1+power((2.*pi.*fdata.*d),2)))+ obj.Sff_sys_noise{kidn,Pindex,nT};
        end
    
       
        function genFknee(obj,kidn,Pindex,nT)
            %genFknee method finds the Fknee based on the fits.
            %p = obj.findp(kidn,Pindex);
            % Analytic 
            
%             obj.Fknee{kidn,Pindex,nT} = power(10,log10(obj.CBfit.Cgr{kidn,Pindex,nT}/obj.TLSfit.Ctls{kidn,Pindex,nT})/(-1*obj.TLSfit.gamma{kidn,Pindex,nT}));
            obj.Fknee{kidn,Pindex,nT} = power(obj.TLSfit.Ctls{kidn,Pindex,nT}/obj.CBfit.Cgr{kidn,Pindex,nT},(1/obj.TLSfit.gamma{kidn,Pindex,nT}));

            obj.Sknee{kidn,Pindex,nT} = obj.CBfit.Cgr{kidn,Pindex,nT}
            
            
            
%             UpA_filled(obj,kidn,Pindex,nT);
%             try
%             [obj.Fknee{kidn,Pindex,nT},obj.Sknee{kidn,Pindex,nT}] = findintersect_SdB2(obj.fTLS_filled,obj.fCB_filled,obj.fguess);
%             catch e
%             fprintf(1,'There was an error! The message was:\n%s',e.message);
%             warning('Problem using function.  Skipping... Setting to NaN');
%             obj.Fknee{kidn,Pindex,nT} = NaN;
%             obj.Sknee{kidn,Pindex,nT} = NaN;
%             figure
%             axes()
%             plot(linspace(obj.fguess(1),obj.fguess(2),1000),obj.fTLS_filled(linspace(obj.fguess(1),obj.fguess(2),1000))-obj.fCB_filled(linspace(obj.fguess(1),obj.fguess(2),1000)))
%             figure
%             axes('XScale','log','YScale','log')
%             hold on
%             plot(linspace(obj.fguess(1),obj.fguess(2),1000),obj.fTLS_filled(linspace(obj.fguess(1),obj.fguess(2),1000))-obj.fCB_filled(linspace(obj.fguess(1),obj.fguess(2),1000)))
%             plot(linspace(obj.fguess(1),obj.fguess(2),1000),-obj.fTLS_filled(linspace(obj.fguess(1),obj.fguess(2),1000))+obj.fCB_filled(linspace(obj.fguess(1),obj.fguess(2),1000)))
%             hold off
%             end
            
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
        
        function plottest(obj,fig_n,ax_n)
            figure(obj.fig(fig_n))
            plot(obj.ax(ax_n),[1 2 3],[4 5 6])
        end
        
        
        function plotsingle(obj,fig_n,ax_n,kidn,Pindex,nT,SW,marker,Colorv)
            %plots to figure currently in focus.
            p = obj.findp(kidn,Pindex);
            UpA(obj,kidn,Pindex,nT)
            obj.Sff{kidn,Pindex,nT} = obj.NOISE(p).FFTnoise{nT}(:,4);
            toplot = obj.Sff{kidn,Pindex,nT} > 0;
            figure(obj.fig(fig_n))
            if SW.plotdata
            plot(obj.ax(ax_n),obj.freq(toplot),lintodb(obj.Sff{kidn,Pindex,nT}(toplot)));
            end
            %plot()
            xlabel('F [Hz]');ylabel('S_F/F^2 [dBc/Hz]')
            %Old: %xlim([0.5,1e5]);grid on;ylim([-220,-140])
            xlim([0.1,5e5]);grid on;ylim([-220,-140])
            title(append('KID#',string(obj.NOISE(p).KIDnumber)," |Power ",string(obj.NOISE(p).ReadPower),"dBm")); 
            %+++++|Begin:Add TLS line--------------------------------------
            if SW.plottls
            %disp(string(size(obj.TLSfit.C)))
            %disp(string(obj.TLSfit.C{kidn,Pindex,nT}(1)))
            plot(obj.ax(ax_n),obj.freq(toplot),lintodb(obj.fTLS(obj.TLSfit.C{kidn,Pindex,nT},obj.freq(toplot))),'-.','Color','black');
            end
            %-------/End:Add TLS line--------------------------------------
            %+++++|Begin:Plot CB fit---------------------------------------
			if SW.plotgr
            plot(obj.ax(ax_n),obj.freq(toplot),lintodb(obj.fCB(obj.CBfit.C{kidn,Pindex,nT},obj.freq(toplot))),':','Color','black','LineWidth',1.5)
            end
            %-------/End:Plot CB fit---------------------------------------
			%+++++|Begin:Visualize Fknee-----------------------------------
			if SW.plotFknee
            plot(obj.ax(ax_n),obj.Fknee{kidn,Pindex,nT},lintodb(obj.Sknee{kidn,Pindex,nT}),marker,'Color',Colorv,'LineWidth',2)
            end
            %-------/End:Visualize Fknee-----------------------------------

            
            
            
            
            %+++++|Begin:fitting ranges------------------------------------
			xline(obj.freq(obj.TLSfit.mini),'--','Color','c','LineWidth',1,'HandleVisibility','off')
            xline(obj.freq(obj.TLSfit.maxi),'--','Color','c','LineWidth',1,'HandleVisibility','off')
            xline(obj.freq(obj.CBfit.mini),'--','Color','m','LineWidth',1,'HandleVisibility','off')
            xline(obj.freq(obj.CBfit.maxi),'--','Color','m','LineWidth',1,'HandleVisibility','off')
            %-------/End:fitting ranges------------------------------------

            
            
        end
        function plotcomptau(obj)
            %
            
        end
        function plotmulti(obj)
            %
            
        end
        %<<<<<End Plotting fuctions.---------------------------------------
        
        %+++++|Begin:Auxillary func----------------------------------------
        function p = findp(obj,kidn,Pindex)% this tells you what p belongs to what power index.
            p = obj.IndexPsort{kidn,1}(Pindex);
        end
        
        
        
		%-------/End:Auxillary func----------------------------------------

        
        
        
        
    end
end

