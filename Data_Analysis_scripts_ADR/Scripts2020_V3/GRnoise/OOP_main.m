%% Refactor OOP version!
clc
clear all
close all

%% (1) Single plots
figurepath = ['..' filesep '..' filesep '..' filesep 'Export_Figures_noGit\OOP\automatisch' filesep ];
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
Single_plot_loop_all = 0;
o = Cfit(FFTsubsubdir,'NOISE_2D.mat');

%-------Do all analysis..
for kidn = kidn_iter
for Pindex= P_iter
for Tindex= T_iter
[CTLS,gamma]= o.genfitsingle(kidn,Pindex,Tindex);
o.genFknee(kidn,Pindex,Tindex);
end
end 
end

%% Temperature variation single plots
o = o.init_figax(fignum,axnum,'loglin');
set(gca, 'FontName', 'Arial')
kidn_iter = 6;
P_iter = 7; 
T_iter = 1:2:14;
Tcolors = colormapJetJB(14);
Pcolors = colormapcoolSdB(7);
%-------
SW.plotdata = 1;
SW.plottls = 1;
SW.plotgr = 1;
SW.plotFknee = 1;
SW.plotTotal = 1;
handleVisible = [{'on'},{'on'},{'on'},{'on'},{'on'}];
for kidn = kidn_iter
for Pindex= P_iter
for Tindex= T_iter
Colorcell = genColorcell(T_iter,Tindex,Tcolors,P_iter,Pindex,Pcolors);
[f1,ax1] = o.plotsingle(fignum,axnum,kidn,Pindex,Tindex,SW,'o',Colorcell,handleVisible);
disp(Pindex)
title(sprintf('$KID_{N}$: %i | P_{read}: %3.0f dBm | $$T_{bath}$$ =%1.3f K | ID:(%i,%i,%i)',kidn,o.getPread(kidn,Pindex,Tindex),o.getT(kidn,Pindex,Tindex),kidn,Pindex,Tindex),'Interpreter','latex')
%legendstr = ['']
end
end 
end
legend('Data','TLS fit','GR + sys noise','Combined','$F_{knee}$','Interpreter','latex')

set(gca,'TickLabelInterpreter', 'latex')
set(gcf,'units','centimeters','position',[20,1,20,0.65*20])
filename = sprintf('KID%iP%3.0fT%1.3f',kidn,o.getPread(kidn,Pindex,Tindex), o.getT(kidn,Pindex,Tindex));
saveas(gca,[figurepath  filename '.fig' ])
exportgraphics(gca,[figurepath  filename '.png' ]) 
exportgraphics(gca,[figurepath filename '.pdf' ])


%% Temperature variation single plots - Loop over all for 1 time - Last run: 14-08 17:52
figurepath = ['..' filesep '..' filesep '..' filesep 'Export_Figures_noGit\OOP\automatisch\single_plot_loop_all_rejection\' filesep ];
if Single_plot_loop_all 
Tcolors = colormapJetJB(14);
Pcolors = colormapcoolSdB(7);

kidn_iter = 1:6;
P_iter = 1:7; 
T_iter = 1:14;
fignum=2;
axnum = fignum;
%-------
SW.plotdata = 1;
SW.plottls = 1;
SW.plotgr = 1;
SW.plotFknee = 1;
SW.plotTotal = 1;
SW.cyanfit = 1;
SW.magentafit = 1;
handleVisible = [{'on'},{'on'},{'on'},{'on'},{'on'}];
for kidn = kidn_iter
for Pindex= P_iter
for Tindex= T_iter
o = o.init_figax(fignum,axnum,'loglin');
figure(o.fig(fignum))
Colorcell = genColorcell(T_iter,Tindex,Tcolors,P_iter,Pindex,Pcolors);
[f1,ax1] = o.plotsingle(fignum,axnum,kidn,Pindex,Tindex,SW,'o',Colorcell,handleVisible);
set(f1,'visible','off')
disp(Pindex)
title(sprintf('$KID_{N}$: %i | $P_{read}$: %3.0f dBm | $$T_{bath}$$ =%1.3f K',kidn,o.getPread(kidn,Pindex,Tindex),o.getT(kidn,Pindex,Tindex)),'Interpreter','latex')
%legendstr = ['']
legend('Data','TLS fit','GR + sys noise','Combined','$F_{knee}$','Interpreter','latex')
set(ax1,'TickLabelInterpreter', 'latex')
set(f1,'units','centimeters','position',[20,1,20,0.65*20])
filename = sprintf('KID%iP%3.0fT%1.3f',kidn,o.getPread(kidn,Pindex,Tindex), o.getT(kidn,Pindex,Tindex));
saveas(gca,[figurepath  filename '.fig' ])
exportgraphics(gca,[figurepath  filename '.png' ]) 
exportgraphics(gca,[figurepath filename '.pdf' ])
close(f1)
end
end 
end
end




%% Temperature variation 4 plot
o = o.init_figax(fignum,axnum,'loglin');
set(gca, 'FontName', 'Arial')
kidn_iter = 1;
P_iter = 7; 
T_iter = [1 6 10 14];
Tcolors = colormapJetJB(14);
Pcolors = colormapcoolSdB(7);
%-------
SW.plotdata = 1;
SW.plottls = 1;
SW.plotgr = 1;
SW.plotFknee = 1;
SW.plotTotal = 1;
SW.cyanfit = 1;
SW.magentafit = 1;
handleVisible = [{'on'},{'off'},{'off'},{'off'},{'off'}];
legendTstr = [];
for kidn = kidn_iter
for Pindex= P_iter
for Tindex= T_iter
Colorcell = genColorcell(T_iter,Tindex,Tcolors,P_iter,Pindex,Pcolors);
[f1,ax1] = o.plotsingle(fignum,axnum,kidn,Pindex,Tindex,SW,'o',Colorcell,handleVisible);
legendTstr = [legendTstr {append(sprintf('%1.3f',o.getT(kidn,Pindex,Tindex)),' K')}];
end
end 
end
legend(legendTstr,'location','eastOutside');
title(sprintf('$KID_{N}$: %i | $P_{read}$: %3.0f dBm',kidn,o.getPread(kidn,Pindex,Tindex)),'Interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex')
set(gcf,'units','centimeters','position',[20,1,20,0.65*20])
filename = '4Temp';
saveas(gca,[figurepath  filename '.fig' ])
exportgraphics(gca,[figurepath  filename '.png' ]) 
exportgraphics(gca,[figurepath filename '.pdf' ])

%% (2) Fknee compare G3,G4 - Temperature
%kidn_iter = 1;
P_iter = 7; 
T_iter = 1:14;
%-------
SW.plotdata = 0;
SW.plottls = 0;
SW.plotgr = 0;
SW.plotFknee = 1;
SW.plotTotal = 0;
SW.cyanfit = 0;
SW.magentafit = 0;
%-------
fignum=2;
axnum = fignum;
Tcolors = colormapJetJB(14);
handleVisible{1} = [{'off'},{'off'},{'off'},{'off'},{'on'}];
handleVisible{2} = [{'off'},{'off'},{'off'},{'off'},{'off'}];
handleVisible{3} = [{'off'},{'off'},{'off'},{'off'},{'off'}];
legendTstr ={};
legendTstr2 ={};
legendTstr3 ={};
o = o.init_figax(fignum,axnum,'loglin');
for kidn = 1:3 % G3
for Pindex= P_iter
for Tindex= T_iter  
Colorcell = genColorcell(T_iter,Tindex,Tcolors,P_iter,Pindex,Pcolors);
[f2,ax2] = o.plotsingle(fignum,axnum,kidn,Pindex,Tindex,SW,'x',Colorcell,handleVisible{kidn});
if kidn==1
legendTstr2 = [legendTstr2 {append(sprintf('G3(2-2-2): %1.3f',o.getT(kidn,Pindex,Tindex)),' K')}];
end
end
end
end

for kidn = 4:6 %G4
for Pindex= P_iter
for Tindex= T_iter
Colorcell = genColorcell(T_iter,Tindex,Tcolors,P_iter,Pindex,Pcolors);
[f2,ax2] = o.plotsingle(fignum,axnum,kidn,Pindex,Tindex,SW,'o',Colorcell,handleVisible{kidn-3});
if kidn==4
legendTstr3 = [legendTstr3 {append(sprintf('G4(4-4-4):%1.3f',o.getT(kidn,Pindex,Tindex)),' K')}];
end
end
end
end
legendTstr = [legendTstr2 legendTstr3];
legend(legendTstr,'location','eastOutside','FontSize',7);
xlim([10,1e4]);grid on;ylim([-190,-160])
title(sprintf('G3 vs G4 | $P_{read}$: Max'),'Interpreter','latex')
filename = 'G3G4compareT';
saveas(gca,[figurepath  filename '.fig' ])
exportgraphics(gca,[figurepath  filename '.png' ]) 
exportgraphics(gca,[figurepath filename '.pdf' ]) 












%% (2) Fknee compare G3,G4 - Power
kidn_iter = 1;
P_iter = 1:7; 
T_iter = 3;
%-------
SW.plotdata = 0;
SW.plottls = 0;
SW.plotgr = 0;
SW.plotFknee = 1;
SW.plotTotal = 0;
SW.cyanfit = 1;
SW.magentafit = 1;
%-------
fignum=3;
axnum = fignum;
Pcolors = colormapcoolSdB(7);
handleVisible = [{'on'},{'on'},{'on'},{'on'},{'on'}];
o = o.init_figax(fignum,axnum,'loglin');
for kidn = 1:3 % G3
for Pindex= P_iter
for Tindex= T_iter
Colorcell = genColorcell(T_iter,Tindex,Tcolors,P_iter,Pindex,Pcolors);
[f3,ax3] = o.plotsingle(fignum,axnum,kidn,Pindex,Tindex,SW,'x',Colorcell,handleVisible);
end
end
end

for kidn = 4:6 %G4
for Pindex= P_iter
for Tindex= T_iter
Colorcell = genColorcell(T_iter,Tindex,Tcolors,P_iter,Pindex,Pcolors);
o.plotsingle(fignum,axnum,kidn,Pindex,Tindex,SW,'o',Colorcell,handleVisible)
end
end
end
filename = 'G3G4compareP';
saveas(gca,[figurepath  filename '.fig' ])
exportgraphics(gca,[figurepath  filename '.png' ]) 
exportgraphics(gca,[figurepath filename '.pdf' ]) 




