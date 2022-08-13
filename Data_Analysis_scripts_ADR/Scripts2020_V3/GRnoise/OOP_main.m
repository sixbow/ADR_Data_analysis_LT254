%% Refactor OOP version!

close all

%% (1) Single plots
FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D'];% This is where the

o = Cfit(FFTsubsubdir,'NOISE_2D.mat')
kidn_iter = 1;
P_iter = 1; 
T_iter = 8;
%-------
fignum=1;
axnum = fignum;
%-------
SW.plottls = 1;
SW.plotgr = 1;
SW.plotFknee = 1;
%-------
for kidn = kidn_iter
for Pindex= P_iter
for Tindex= T_iter
[CTLS,gamma]= o.genfitsingle(kidn,Pindex,Tindex);
o.genFknee(kidn,Pindex,Tindex)
end
end 
end
o = o.init_figax(fignum,axnum,'loglin');
for kidn = kidn_iter
for Pindex= P_iter
for Tindex= T_iter
o.plotsingle(fignum,axnum,kidn,pindex,Tindex,SW)
end
end 
end





