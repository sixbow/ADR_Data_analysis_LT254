classdef Cfit %< handle
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
        %Looping vars
        kidn_iter
        power_iter
        Tbath_iter
        %Fitting vars
        TLSf
        CBf
        fguess
        % Anonymous functions
        fTLS
        fCB
        fline
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
            obj.TLSf.min = 0.8;% [Hz]
            obj.TLSf.max= 30;% [Hz]
            [~,obj.TLSf.mini] = min(abs(obj.freq-obj.TLSf.min)); % This gives the min index
            [~,obj.TLSf.maxi] = min(abs(obj.freq-obj.TLSf.max)); % This gives the max index
            obj.TLSf.C0 = [(10^-15) 1.5];% [CTLS gamma]  linear scale. CTLS at 1 Hz , power with what  the power law decreases. gamma TLS ~0.5
 			%-------/End:TLS fitting parameters (Default!)-----------------
            %+++++|Begin:Combined fitting parameters (Default!)------------
            obj.CBf.min = 1000;% [Hz]
            obj.CBf.max = 150000;%  [Hz]
            [~,obj.CBf.mini] = min(abs(obj.freq-obj.CBf.min)); % "
            [~,obj.CBf.maxi] = min(abs(obj.freq-obj.CBf.max)); % "
            obj.CBf.C0 = [(10^-18) 1.591549430918954e-06];% CTLS at 1 Hz , power with what  the power law decreases. TLS ~0.5
            obj.CBf.lb = [10^-20 (1/(2*pi*obj.CBf.max))];
            obj.CBf.ub = [10^-15 (1/(2*pi*obj.CBf.min))];
            obj.fguess = 2000; % [Hz] initial guess for the intersection point of the TLS and GR noise. 
			%-------/End:Combined fitting parameters (Default!)------------
            obj.test = 3.14;
            disp('Initialized Cfit object!')
        end % End constructor
        %+++++:Begin fitting and calcuations.-------------------------------
        
        
        function  [Ctls,gamma]= genfitsingle(obj,kidn,pindex,nT)
            %genfit method fits the TLS fit and the GR fits 
            p = obj.findp(kidn,pindex);
            obj.Sff{p,nT} = obj.NOISE(p).FFTnoise{nT}(:,4);
            setup_noise_eval_freq = 350000;% [Hz]
            [~,Snoise_i] = min(abs(obj.freq-setup_noise_eval_freq)); % "
            obj.Sff_sys_noise{p,nT} = obj.NOISE(p).FFTnoise{nT}(Snoise_i,4); % This is the setup noise level evaluated at 350000Hz
            %Updating and Setting anonymous functions
            obj.UpA(p,nT)
			%+++++|Begin:TLS fitting---------------------------------------
            dBlogDataX = log10(obj.freq);
            dBlogDataY = 10.*log10(obj.NOISE(p).FFTnoise{nT}(:,4));
            TLS_coof_dBlog{p,nT} = polyfit(dBlogDataX(obj.TLSf.mini:obj.TLSf.maxi),dBlogDataY(obj.TLSf.mini:obj.TLSf.maxi),1);
            obj.TLSf.Ctls{p,nT} =  power(10,(TLS_coof_dBlog{p,nT}(2)/10));
            obj.TLSf.gamma{p,nT} = -(TLS_coof_dBlog{p,nT}(1)/10);
            obj.TLSf.C{p,nT} = [obj.TLSf.Ctls obj.TLSf.gamma];
            Ctls = obj.TLSf.Ctls;% Output to user
            gamma = obj.TLSf.gamma;% Output to user
			%-------/End:TLS fitting---------------------------------------
         
			%+++++|Begin:Combined fitting----------------------------------
			%-------/End:Combined fitting---------------------------------- 
        end% End genfitsingle
        
        function UpA(obj,p,nT)
            obj.fTLS = @(C_v,fdata)C_v(1).*power(fdata,-1*C_v(2)); %Model we use to fit C_v is the constants vector that we are trying to find.
            obj.fCB  = @(C_v,fdata) C_v(1)./((1+power((2.*pi.*fdata.*C_v(2)),2).*(1)))+ obj.Sff_sys_noise{p,nT};
            obj.fline = @(C_v,xdata)C_v(1).*xdata+C_v(2);
        end
    
        
        
        function genFknee(obj)
            %genFknee method finds the Fknee based on the fits.
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
        
        
        function plotsingle(obj,fig_n,ax_n,kidn,pindex,nT)
            %plots to figure currently in focus.
            p = obj.findp(kidn,pindex);
            UpA(obj,p,nT)
            obj.Sff{p,nT} = obj.NOISE(p).FFTnoise{nT}(:,4);
            toplot = obj.Sff{p,nT} > 0;
            figure(obj.fig(fig_n))
            plot(obj.ax(ax_n),obj.freq(toplot),lintodb(obj.Sff{p,nT}(toplot)));
            %plot()
            xlabel('F [Hz]');ylabel('S_F/F^2 [dBc/Hz]')
            %Old: %xlim([0.5,1e5]);grid on;ylim([-220,-140])
            xlim([0.1,5e5]);grid on;ylim([-220,-140])
            title(append('KID#',string(obj.NOISE(p).KIDnumber)," |Power ",string(obj.NOISE(p).ReadPower),"dBm")); 
            
        end
        function plotcomptau(obj)
            %
            
        end
        function plotmulti(obj)
            %
            
        end
        %<<<<<End Plotting fuctions.---------------------------------------
        
        %+++++|Begin:Auxillary func----------------------------------------
        function p = findp(obj,kidn,power_selector)% this tells you what p belongs to what power index.
            p = obj.IndexPsort{kidn,1}(power_selector);
        end
        
        
        
		%-------/End:Auxillary func----------------------------------------

        
        
        
        
    end
end

