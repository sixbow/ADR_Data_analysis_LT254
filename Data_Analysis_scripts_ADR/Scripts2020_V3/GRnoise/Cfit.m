classdef Cfit
    %Cfit is the base class for our data analysis.
    %
    
    properties
        %Data vars
        NOISE
        IndexPsort
        IndexPopt
        freq
        %Looping vars
        kidn_iter
        power_iter
        Tbath_iter
        %Fitting vars
        TLSf
        CBf
        fguess
        %Plotting vars
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
            
            disp('Initialized Cfit object!')
        end % End constructor
        %+++++:Begin fitting and calcuations.-------------------------------
        
        
        
        function genfit(self)
            %genfit method fits the TLS fit and the GR fits 
        end
        
        
        
        
        
        
        
        function genFknee(self)
            %genFknee method finds the Fknee based on the fits.
        end
        %-----:End fitting and calcuations.---------------------------------
        %>>>>>:Begin Plotting fuctions.-------------------------------------
        function plotsingle(self)
            %   
        end
        function plotcomptau(self)
            %
            
        end
        function plotmulti(self)
            %
            
        end
        %<<<<<End Plotting fuctions.---------------------------------------
    end
end

