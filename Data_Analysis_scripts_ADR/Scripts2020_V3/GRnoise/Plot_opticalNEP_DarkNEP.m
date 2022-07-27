function Plot_opticalNEP_DarkNEP(opticalnepisthere)

opticalnepisthere = 0; % default 1, 0 if no optical NEP

addpath([pwd,filesep,'..',filesep,'subroutines']);

if opticalnepisthere == 1
    addpath([pwd,filesep,'..',filesep,'subroutines']);
    OptNEPFile = 'KIDparam.mat';
    OptNEPdir = '/Volumes/kid/KIDonSun/experiments/Entropy ADR/LT179/LT179-chip4-optical/Run 2 - Federica/Combined/2D_BB/2D_BB';
    load([OptNEPdir filesep OptNEPFile],'OptNEPmin');
    for kidn = 1 : length(OptNEPmin.KIDID)
        OptNEPmin.AluLength(kidn) = NEPmin.AluLength(OptNEPmin.KIDID(kidn) == NEPmin.KIDID);
    end
end

DarkNEPFile = 'NEP_P.mat';
DarkNEPdir = '\\MARS\kid\KIDonSun\experiments\Entropy ADR\LT165_W2_Chip10\NEP';
load([DarkNEPdir filesep DarkNEPFile],'NEPmin');

%Calcualte the GR NEP. Volume is already in um^3
for kidn = 1 : length(NEPmin.KIDID)
    [NEPmin.NEPGR(kidn),NEPmin.nqp(kidn),NEPmin.Nqp(kidn)] = NEPquick2(NEPmin.tau(kidn),NEPmin.Volume(kidn));
end


figure(1)
if opticalnepisthere == 1
    %optical
    subplot(2,2,1)
    semilogy(OptNEPmin.KIDID,OptNEPmin.theta,'or','MarkerSize',8,'MarkerfaceColor','r');  hold on
    xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-17]);
    title('Phase readout')
    subplot(2,2,2)
    semilogy(OptNEPmin.KIDID,OptNEPmin.R,'or','MarkerSize',8,'MarkerfaceColor','r');  hold on
    xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-17]);
    title('R readout')
    subplot(2,2,3)
    semilogy(OptNEPmin.AluLength,OptNEPmin.theta,'or','MarkerSize',8,'MarkerfaceColor','r');  hold on
    xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-17]);
    title('Phase readout')
    subplot(2,2,4)
    semilogy(OptNEPmin.AluLength,OptNEPmin.R,'or','MarkerSize',8,'MarkerfaceColor','r');  hold on
    xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-17]);
    title('R readout')
end

%dark
subplot(2,2,1)
semilogy(NEPmin.KIDID,NEPmin.theta,'ob','MarkerSize',8,'MarkerfaceColor','b');  hold on
xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-17]);
title('Phase readout')
subplot(2,2,2)
semilogy(NEPmin.KIDID,NEPmin.R,'ob','MarkerSize',8,'MarkerfaceColor','b');  hold on
xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-17]);
title('R readout')
subplot(2,2,3)
semilogy(NEPmin.AluLength,NEPmin.theta,'ob','MarkerSize',8,'MarkerfaceColor','b');  hold on
xlabel('Alu Length  (\mum)');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-17]);
title('Phase readout')
subplot(2,2,4)
semilogy(NEPmin.AluLength,NEPmin.R,'ob','MarkerSize',8,'MarkerfaceColor','b');  hold on
xlabel('Alu Length  (\mum)');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-17]);
title('R readout')

%GR
subplot(2,2,1)
semilogy(NEPmin.KIDID,NEPmin.NEPGR,'sk','MarkerSize',8,'MarkerfaceColor','k');  hold on
xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-17]);
title('Phase readout')
subplot(2,2,2)
semilogy(NEPmin.KIDID,NEPmin.NEPGR,'sk','MarkerSize',8,'MarkerfaceColor','k');  hold on
xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-17]);
title('R readout')
subplot(2,2,3)
semilogy(NEPmin.AluLength,NEPmin.NEPGR,'sk','MarkerSize',8,'MarkerfaceColor','k');  hold on
xlabel('Alu Length  (\mum)');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-17]);
title('Phase readout')
subplot(2,2,4)
semilogy(NEPmin.AluLength,NEPmin.NEPGR,'sk','MarkerSize',8,'MarkerfaceColor','k');  hold on
xlabel('Alu Length  (\mum)');ylabel('NEP (W /\surd Hz)');grid on;ylim([1e-20 1e-17]);
title('R readout')

for nm = 1:4
    subplot(2,2,nm)
    if opticalnepisthere == 1
        legend('optical NEP','Dark NEP','GR' )
    else
        legend('Dark NEP','GR' )
    end
    
    MakeGoodFigure(12,11,12,[DarkNEPdir filesep  'NEPcompare'])
end
end

