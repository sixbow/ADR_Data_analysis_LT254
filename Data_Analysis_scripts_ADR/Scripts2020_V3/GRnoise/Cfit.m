classdef Cfit
    %Cfit is the base class for our data analysis.
    %
    
    properties
        NOISE
        IndexPsort
        IndexPopt
        freq
        kidn_iter
        power_iter
        Tbath_iter
        %SW
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
            disp('Initialized Cfit object!')
            obj.freq = obj.NOISE(1).FFTnoise{1,1}(:,1);
            obj.kidn_iter = 1:6; % Default value
            obj.power_iter = obj.IndexPsort{obj.kidn_iter,1};
            obj.Tbath_iter = 1:14;
            %obj.SW = SW; % struct that contains switches.
        end
        
        function genfit(self)
            %genfit method fits the TLS fit and the GR fits 
        end
        function genFknee(self)
            %genFknee method finds the Fknee based on the fits.
        end
        function plotsingle(self)
            %
            
        end
        function plotcomptau(self)
            %
            
        end
        function plotmulti(self)
            %
            
        end
    end
end

