% This code is used as a subroutine that changes
% CrossPSDNoise.mat to remove the offset of the TLS noise in the signal to
% improve fitting
% TODO: Rename  CrossPSDNoise_2D -> CrossPSDNoise_2D_INPUT
% Step 1: importing the mat file 
% Step 2: getting the offset of the TLS-noise
% Step 3: Subtract and save
% Step 4: Backup old imported file 
% Author: Sietse de Boer 
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





save([Outputfolderdir,filesep,matfile],'NOISE','IndexP_sub_opt','KIDnumbers');
save([Outputfolderdir,filesep,matfile2],'CrossPSDNOISE');
disp('Jobs done!');
