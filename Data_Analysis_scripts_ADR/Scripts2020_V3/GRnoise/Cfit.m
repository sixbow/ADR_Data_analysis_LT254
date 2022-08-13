classdef Cfit
    %Cfit is the base class for our data analysis.
    %
    
    properties
        NOISE
        IndexPsort
        IndexPopt
    end
    
    methods
        function obj = Cfit(FFTsubsubdir,filename)
            %Cfit Construct an instance of this class
            addpath('..\subroutines')
            ChipInfo_path = ['..' filesep '..' filesep ]; %root path where data is, one higher than the scripts       
            obj.NOISE = load([ChipInfo_path FFTsubsubdir filesep filename],'NOISE').NOISE;
            obj.IndexPsort = load([ChipInfo_path FFTsubsubdir filesep filename],'IndexPsort').IndexPsort;
            obj.IndexPopt = load([ChipInfo_path FFTsubsubdir filesep filename],'IndexPopt').IndexPopt;
            disp('Initialized Cfit object!')
            disp(size(obj.IndexPopt));
        end
        
        function outputArg = gettest()
            %METHOD1 Summary of this method goes here
            disp(size(obj.IndexPopt));
        end
    end
end

