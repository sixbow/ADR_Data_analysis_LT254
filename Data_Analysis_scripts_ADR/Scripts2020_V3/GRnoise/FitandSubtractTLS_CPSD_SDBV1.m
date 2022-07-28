% This code is used as a subroutine that changes
% CrossPSDNoise.mat to remove the offset of the TLS noise in the signal to
% improve fitting
% TODO: Rename  CrossPSDNoise_2D -> CrossPSDNoise_2D_INPUT
% Step 1: importing the mat file 
% Step 2: getting the offset of the TLS-noise
% Step 3: Subtract and save
% Step 4: Backup old imported file 
% Author: Sietse de Boer 
clc
clear all;%This should be commented out for performance reasons.
close all;
ChipInfo_path = ['..' filesep '..' ]; %root path where data is, one higher than the scripts
FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D'];% This is where the
Outputfolder = 'CPSDMinusTLS';
Outputfolderdir = [ChipInfo_path,filesep,FFTsubsubdir,filesep,Outputfolder] ; 
% CrossPSDNoise_2D file is. %FFTsubdir = [filesep 'Noise_Powers_165mK' filesep 'FFT' filesep 'Power'];     %
% Loading files!
matfile = 'Noise_2D.mat';
matfile2 = 'CrossPSDNoise_2D';
matfile3 = 'CrossPSDFit_2D';
load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','IndexP_sub_opt','KIDnumbers');
load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile2],'CrossPSDNOISE');
CPSDfolder = [ChipInfo_path, filesep, 'CrossPSD2D'] ;
if ~exist(Outputfolderdir, 'dir')
       mkdir(Outputfolderdir)
end % This makes a subfolder called CPSDMinusTLS which will be were the files for this script will be created
% will be making for loop, But for now 
%%
%f1 = figure();
%hold on
kidn = 1;
nT = 1;%Likely loop between 1:14  ~length()
p = 1; % need to loop over all powers.
%Data
Current_freq = CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1);
Current_S_CPSD = abs(real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2)));
%ax1 = axes('XScale', 'log', 'YScale', 'log');
loglog(Current_freq,Current_S_CPSD);
% x(1) = C_{TLS}
% x(2) = \alpha ~0.5
% model S_{TLS} = C_{TLS}*(f)^{\alpha}
%begin_data_point = 1; 
%end_data_point = 30; % We only want to fit the first part of the TLS noise since there it is dominant
% we have about 50 points so first 5 points seems fine
%x0 = [(10^-7) 0.5];
%Model_TLS = @(x,fdata)x(1)*(fdata).^(x(2));

%[x,resnorm,~,exitflag,output] = lsqcurvefit(Model_TLS,x0,Current_freq(begin_data_point:end_data_point),Current_S_CPSD(begin_data_point:end_data_point))
%loglog(Current_freq,Model_TLS(x,Current_freq));
% xline(Current_freq(begin_data_point),'--')
% xline(Current_freq(end_data_point),'--')
% f2 = figure;
% plot(Current_freq,Current_S_CPSD)
% plot(Current_freq,Model_TLS(x,Current_freq));
% xline(Current_freq(begin_data_point),'--')
% xline(Current_freq(end_data_point),'--')

% \will be making for Power_loop, But for now <--- Place end here!
% \will be making for Temp_loop, But for now <--- Place end here!
%hold off
% \will be making for KID_loop, But for now <--- Place end here!

save([Outputfolderdir,filesep,matfile],'NOISE','IndexP_sub_opt','KIDnumbers');
save([Outputfolderdir,filesep,matfile2],'CrossPSDNOISE');
disp('Jobs done!');
% disp('Press any key to quit and clean up!')
% pause
% clear all