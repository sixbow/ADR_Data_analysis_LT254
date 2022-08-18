%% Refactor OOP version!
clc
clear all
close all

%%  Paths..
figurepath = ['..' filesep '..' filesep '..' filesep 'Export_Figures_noGit\OOP\automatisch' filesep ];
FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D'];% This is where the
addpath('..\subroutines')
%% (1) Analysis
SW_simsome = 0;% Ability to choose some values.(Default:0)
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

SW.CrossTau =1 ; %1 Uses JBfitted Crosstau_qp, 0 Fits from frac freq PSD.
o = Cfit(FFTsubsubdir,'NOISE_2D.mat','CrossPSDFit_2D.mat');
o.genCross_tau(kidn_iter,P_iter,T_iter)  % generates the nessesary entries for the cross PSD fitted tau. 
%-------Do all analysis..
for kidn = kidn_iter
for Pindex= P_iter
for Tindex= T_iter
if SW.CrossTau 
[CTLS,gamma]= o.genfitCrosssingle(kidn,Pindex,Tindex);
else
[CTLS,gamma]= o.genfitsingle(kidn,Pindex,Tindex);
end
o.genfRing(kidn,Pindex,Tindex)
o.genFknee(kidn,Pindex,Tindex);
end
end 
end


%% Temperature variation single plots
o = o.init_figax(fignum,axnum,'loglin');
SW.usePopt = 0;% Use the 120mK Popt value. 
set(gca, 'FontName', 'Arial')
kidn_iter = 1;
% Write a fuction that creates Pindex = o.findPopt(kidn)
P_iter = o.findPiopt(kidn_iter); 
T_iter = T_iter+1;
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
SW.plotringline = 1;
handleVisible = [{'on'},{'on'},{'on'},{'on'},{'on'},{'on'}];
for kidn = kidn_iter
for Pindex= P_iter
for Tindex= T_iter
    
Colorcell = genColorcell(T_iter,Tindex,Tcolors,P_iter,Pindex,Pcolors);
[f1,ax1] = o.plotsingle(fignum,axnum,kidn,Pindex,Tindex,SW,'o',Colorcell,handleVisible);
disp(Pindex)
title(sprintf('$KID_{N}$: %i | $P_{read}$: %3.0f dBm | $$T_{bath}$$ =%1.3f K | ID:(%i,%i,%i)',kidn,o.getPread(kidn,Pindex,Tindex),o.getT(kidn,Pindex,Tindex),kidn,Pindex,Tindex),'Interpreter','latex')
%legendstr = ['']
end
end 
end
legend('Data','TLS fit','GR + sys noise','Combined','$F_{knee}$','$f_{ring}$','Interpreter','latex')

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
handleVisible = [{'on'},{'on'},{'on'},{'on'},{'on'},{'off'}];
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
handleVisible = [{'on'},{'off'},{'off'},{'off'},{'off'},{'off'}];
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
handleVisible{1} = [{'off'},{'off'},{'off'},{'off'},{'on'},{'off'}];
handleVisible{2} = [{'off'},{'off'},{'off'},{'off'},{'off'},{'off'}];
handleVisible{3} = [{'off'},{'off'},{'off'},{'off'},{'off'},{'off'}];
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


%% Fknee vs T

Pindex = 7;
%--------
fignum=2;
axnum = fignum;
o = o.init_figax(fignum,axnum,'linlog');
for kidn=1:6
o.plottempfknee(fignum,axnum,kidn,Pindex)
end
grid on
legend('KID:1','KID:2','KID:3','KID:4','KID:5','KID:6')
title('$T_{bath}$ vs $f_{knee} $@$P_{read}=$Max','Interpreter','latex')
xlabel('$T_{bath}$ [K]','Interpreter','latex')
ylabel('$f_{knee}$ [Hz]','Interpreter','latex')
filename = sprintf('FkneevsTKID%i',kidn);
saveas(gca,[figurepath  filename '.fig' ])
exportgraphics(gca,[figurepath  filename '.png' ]) 
exportgraphics(gca,[figurepath filename '.pdf' ]) 





%% Fknee vs Qi


Pindex = 7;
%--------
fignum=3;
axnum = fignum;
o = o.init_figax(fignum,axnum,'loglog');
for kidn=1:6
o.plotQifknee(fignum,axnum,kidn,Pindex)
end
grid on
legend('KID:1','KID:2','KID:3','KID:4','KID:5','KID:6')
title('$Q_{i}$ vs $f_{knee}$ @$P_{read}=$Max','Interpreter','latex')
xlabel('$Q_{i}$','Interpreter','latex')
ylabel('$f_{knee}$ [Hz]','Interpreter','latex')
filename = sprintf('FkneevsQiKID%i',kidn);
saveas(gca,[figurepath  filename '.fig' ])
exportgraphics(gca,[figurepath  filename '.png' ]) 
exportgraphics(gca,[figurepath filename '.pdf' ]) 

%% Tau vs Temp
Tcolors = colormapJetJB(14);
for kidn = 1:6
Pindex = 7
T_iter = 1:14
%kidn = 6
fignum=4;
axnum = fignum;
o = o.init_figax(fignum,axnum,'linlog');

o.plottautemp(fignum,axnum,kidn,Pindex,Tcolors,6)
grid on
title(sprintf('$KID_{N}$ %i | T vs $\\tau_{qp}$ @ $P_{read} =$ Max',kidn),'Interpreter','latex')
xlabel('$T_{bath} [K]$','Interpreter','latex')
ylabel('$\tau_{qp}$ [msec]','Interpreter','latex')

filename = sprintf('TauvsTempKID%i',kidn);
saveas(gca,[figurepath  filename '.fig' ])
exportgraphics(gca,[figurepath  filename '.png' ]) 
exportgraphics(gca,[figurepath filename '.pdf' ])
end


%% (2) Fknee compare G3,G4 - Power
kidn_iter = 1;
P_iter = 1:7; 
T_iter = 13;
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
handleVisible = [{'on'},{'on'},{'on'},{'on'},{'on'},{'off'}];
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

%% Power variation 4 plot
o = o.init_figax(fignum,axnum,'loglin');
set(gca, 'FontName', 'Arial')
kidn_iter = 1;
P_iter = 1:7; 
T_iter = 13;
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
handleVisible = [{'on'},{'off'},{'off'},{'off'},{'off'},{'off'}];
legendTstr = [];
for kidn = kidn_iter
for Pindex= P_iter
for Tindex= T_iter
Colorcell = genColorcell(T_iter,Tindex,Tcolors,P_iter,Pindex,Pcolors);
[f1,ax1] = o.plotsingle(fignum,axnum,kidn,Pindex,Tindex,SW,'o',Colorcell,handleVisible);
legendTstr = [legendTstr {append(sprintf('%1.3f',o.getPread(kidn,Pindex,Tindex)),' dBm')}];
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



