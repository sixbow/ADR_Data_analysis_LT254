%% Refactor OOP version!

close all

%% (1) Single plots
figurepath = ['..' filesep '..' filesep '..' filesep 'Export_Figures_noGit\OOP\automatisch' filesep ]
FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D'];% This is where the
addpath('..\subroutines')
%% Plotting: single plot!
SW_simsome = 0;
if SW_simsome
kidn_iter = 1;
P_iter = 6; 
T_iter = 12;
elseif SW_simsome == 0
kidn_iter = 1:6;
P_iter = 1:7; 
T_iter = 1:14;
end
%-------
fignum=1;
axnum = fignum;
%-------
SW.plotdata = 0;
SW.plottls = 0;
SW.plotgr = 0;
SW.plotFknee = 1;

o = Cfit(FFTsubsubdir,'NOISE_2D.mat')

%-------Do all analysis..
for kidn = kidn_iter
for Pindex= P_iter
for Tindex= T_iter
[CTLS,gamma]= o.genfitsingle(kidn,Pindex,Tindex);
o.genFknee(kidn,Pindex,Tindex)
end
end 
end

%% Temperature variation single plots
o = o.init_figax(fignum,axnum,'loglin');
kidn_iter = 1;
P_iter = 7; 
T_iter = 13;
%-------
SW.plotdata = 1;
SW.plottls = 1;
SW.plotgr = 1;
SW.plotFknee = 1;
for kidn = kidn_iter
for Pindex= P_iter
for Tindex= T_iter
o.plotsingle(fignum,axnum,kidn,Pindex,Tindex,SW,'o','black')
end
end 
end
legend('Data','TLS fit','GR + sys noise','Combined')




%% Fknee compare G3,G4
P_iter = 7; 
T_iter = 1:14;
%-------
SW.plotdata = 0;
SW.plottls = 0;
SW.plotgr = 0;
SW.plotFknee = 1;
%-------
fignum=2;
axnum = fignum;
Tcolors = colormapJetJB(14);
o = o.init_figax(fignum,axnum,'loglin');
for kidn = 1:3 % G3
for Pindex= P_iter
for Tindex= T_iter
o.plotsingle(fignum,axnum,kidn,Pindex,Tindex,SW,'x',Tcolors(Tindex,:))
end
end
end

for kidn = 4:6 %G4
for Pindex= P_iter
for Tindex= T_iter
o.plotsingle(fignum,axnum,kidn,Pindex,Tindex,SW,'o',Tcolors(Tindex,:))
end
end
end
filename = 'G3G4compare';
saveas(gca,[figurepath  filename '.fig' ])
exportgraphics(gca,[figurepath  filename '.png' ]) 
exportgraphics(gca,[figurepath filename '.pdf' ]) 






%% Fknee compare G3,G4
P_iter = 1:7; 
T_iter = 3;
%-------
SW.plotdata = 0;
SW.plottls = 0;
SW.plotgr = 0;
SW.plotFknee = 1;
%-------
fignum=2;
axnum = fignum;
Pcolors = colormapcoolSdB(7);
o = o.init_figax(fignum,axnum,'loglin');
for kidn = 1:3 % G3
for Pindex= P_iter
for Tindex= T_iter
o.plotsingle(fignum,axnum,kidn,Pindex,Tindex,SW,'x',Pcolors(Pindex,:))
end
end
end

for kidn = 4:6 %G4
for Pindex= P_iter
for Tindex= T_iter
o.plotsingle(fignum,axnum,kidn,Pindex,Tindex,SW,'o',Pcolors(Pindex,:))
end
end
end
filename = 'G3G4compareP';
saveas(gca,[figurepath  filename '.fig' ])
exportgraphics(gca,[figurepath  filename '.png' ]) 
exportgraphics(gca,[figurepath filename '.pdf' ]) 




