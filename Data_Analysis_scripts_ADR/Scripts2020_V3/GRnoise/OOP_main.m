%% Refactor OOP version!

close all

%% (1) Single plots
FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D'];% This is where the

o = Cfit(FFTsubsubdir,'NOISE_2D.mat')
kidn = 1;
pindex = 7; 
Tindex = 10;
%-------
fignum=1;
axnum = fignum;
%-------
SW.plottls = 1;
SW.plotgr = 1;
SW.plotFknee = 1;
%-------
[CTLS,gamma]= o.genfitsingle(kidn,pindex,Tindex);
o.genFknee(kidn,pindex,Tindex)
o = o.init_figax(fignum,axnum,'loglin');
o.plotsingle(fignum,axnum,kidn,pindex,Tindex,SW)







