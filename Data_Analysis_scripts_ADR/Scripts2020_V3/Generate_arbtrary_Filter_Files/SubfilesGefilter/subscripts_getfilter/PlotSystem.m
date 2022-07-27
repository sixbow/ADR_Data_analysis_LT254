function [S21,Power] = PlotSystem(Fr,Area,Filter,name,Fi,Frange,anglein)
%calcualtes powers and plot system
%S21 is total transmission
%power is struct with powers coupled and entering cryostat
if nargin == 6
    anglein = 60;
end
kleur = colormapJetJB(length(Fi));

ddf = Fr(2) - Fr(1);

[Power.tot,NEP]=BBradiator_mulitmode(300, Fr(1)/1e12, Fr(end)/1e12 ,1, Area, anglein);
Irad = interp1(NEP.F,NEP.filterirad,Fr);%irradiation over entire Fr range interpolated int Fr matrix
S21 = totalS21(Filter(Fi));     %S21 of filter stack
S21_nowindow = totalS21(Filter(Fi(2:end)));     %S21 of filter stack
IradS21 = Irad.* S21;       %filtered Irad
IradS21_nowindow = Irad.* S21_nowindow;       %filtered Irad
Power.IradS21 = IradS21;
Power.IradS21_nowindow = IradS21_nowindow;
Power.Irad = Irad; 

Power.Pfullband = sum(Power.IradS21) * ddf; %total power transmitted to cold
Power.Pfullband_nowindow = sum(Power.IradS21_nowindow) * ddf; %total power transmitted to cold
Power.Penter = sum(Power.Irad) * ddf;   %total power entering 
Power.Pinband = sum(Power.IradS21(Frange)) * ddf; %total power transmitted in band
%
subplot(3,1,1)
for n = 1 : length(Fi)
     plot(Fr/1e12,Filter{Fi(n)}(1,:),'Linewidth',2,'color',kleur(n,:));hold on
end
plot(Fr/1e12, S21,'k','linewidth',2);
xlim([0 1]);ylim([0 1]);grid on
ylabel('Transmission');xlabel('Frequency  (THz)')
legend({name{Fi}, 'all'},'Location','SouthEast')

subplot(3,1,2)
for n = 1 : length(Fi)
     loglog(Fr/1e12,Filter{Fi(n)}(1,:),'Linewidth',2,'color',kleur(n,:));hold on
end
loglog(Fr/1e12, S21,'k','linewidth',2);
legend({name{Fi}, 'all'},'Location','SouthWest')
xlim([0.1 100]);ylim([1e-8 1]);grid on
ylabel('Transmission');xlabel('Frequency  (THz)')
title(['R = ' num2str((Area/pi).^0.5*1e3,'%.0f') ' mm, angle = ' num2str(anglein,'%.0f') ' degree']);

subplot(3,1,3)
%loglog(Fr, IradS21,'g','linewidth',2);hold on;%_nowindow
loglog(Fr, IradS21_nowindow,'r','linewidth',2);hold on;%_nowindow
loglog(Fr, Irad,'b','linewidth',2);
%plot(Fr(Frange),Power.IradS21(Frange),'-ok')
%legend('Filtered Irradiance','Filtered Irradiance, exl. window','Planck 300K','Location','SouthWest')
legend('Filtered Irradiance','Planck 300K','Location','SouthWest')
grid on;xlim([100e9 Fr(end)]);ylim([1e-23 1e-15])
title(['In band = ' num2str(100*Power.Pinband/Power.Pfullband,'%.0f') ' %, P_{in} = ' num2str(Power.Penter *1e3,'%.0f')...
    'mW, P_T = ' num2str(Power.Pfullband_nowindow *1e3,'%.3f') ' mW' ])
ylabel('Spectral power  (W/Hz)');xlabel('Frequency  (Hz)')
end