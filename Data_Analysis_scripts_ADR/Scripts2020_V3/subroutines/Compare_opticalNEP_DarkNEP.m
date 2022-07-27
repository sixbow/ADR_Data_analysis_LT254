function Compare_opticalNEP_DarkNEP(NEPmin,OptNEPdir)

OptNEPFile = 'KIDparam.mat';
if nargin == 2 && ~isempty(OptNEPdir)
    opticalnepisthere = 1;
    load([OptNEPdir filesep OptNEPFile],'OptNEPmin');
    for kidn = 1 : length(OptNEPmin.KIDID)
        indtouse = OptNEPmin.KIDID(kidn) == NEPmin.KIDID;
        if sum(indtouse) ~= 0
            OptNEPmin.AluLength(kidn) = NEPmin.AluLength(OptNEPmin.KIDID(kidn) == NEPmin.KIDID);
        else
            OptNEPmin.AluLength(kidn) = NaN;
        end
    end
else 
    opticalnepisthere = 0;
    OptNEPdir = '';
end


figure(1000)
if opticalnepisthere == 1
    %optical
    subplot(2,3,1)
    semilogy(OptNEPmin.KIDID,OptNEPmin.theta_Pabs,'or','MarkerSize',8,'MarkerfaceColor','r');  hold on
    subplot(2,3,2)
    semilogy(OptNEPmin.KIDID,OptNEPmin.ddxdPabs,'or','MarkerSize',8,'MarkerfaceColor','r');  hold on
    subplot(2,3,3)
    semilogy(OptNEPmin.KIDID,OptNEPmin.R_Pabs,'or','MarkerSize',8,'MarkerfaceColor','r');  hold on
    %vs AluLength
    subplot(2,3,4)
    semilogy(OptNEPmin.AluLength,OptNEPmin.theta_Pabs,'or','MarkerSize',8,'MarkerfaceColor','r');  hold on
    subplot(2,3,5)
    semilogy(OptNEPmin.AluLength,OptNEPmin.ddxdPabs,'or','MarkerSize',8,'MarkerfaceColor','r');  hold on
        subplot(2,3,6)
    semilogy(OptNEPmin.AluLength,OptNEPmin.R_Pabs,'or','MarkerSize',8,'MarkerfaceColor','r');  hold on

end

%dark NEPmin.dxdPdark
subplot(2,3,1)
semilogy(NEPmin.KIDID,NEPmin.theta,'ob','MarkerSize',6,'MarkerfaceColor','b');  hold on
xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
title('Phase readout')
subplot(2,3,2)
semilogy(NEPmin.KIDID,NEPmin.dxdPdark,'ok','MarkerSize',8,'MarkerfaceColor','k');  hold on
semilogy(NEPmin.KIDID,1e20*NEPmin.dxdN,'ob','MarkerSize',8,'MarkerfaceColor','b');  hold on
xlabel('KID ID');ylabel('dx/dP  (1/W)');grid on;
subplot(2,3,3)
semilogy(NEPmin.KIDID,NEPmin.R,'ob','MarkerSize',6,'MarkerfaceColor','b');  hold on
xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
title('R readout')

%vs Alu length
subplot(2,3,4)
semilogy(NEPmin.AluLength,NEPmin.theta,'ob','MarkerSize',6,'MarkerfaceColor','b');  hold on
xlabel('Alu Length  (\mum)');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
title('Phase readout')
subplot(2,3,5)
semilogy(NEPmin.AluLength,NEPmin.dxdPdark,'ok','MarkerSize',8,'MarkerfaceColor','k');  hold on
semilogy(NEPmin.AluLength,1e20*NEPmin.dxdN,'ob','MarkerSize',8,'MarkerfaceColor','b');  hold on
xlabel('Alu Length  (\mum)');ylabel('dx/dP  (1/W)');grid on;
subplot(2,3,6)
semilogy(NEPmin.AluLength,NEPmin.R,'ob','MarkerSize',6,'MarkerfaceColor','b');  hold on
xlabel('Alu Length  (\mum)');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
title('R readout')

%GR
subplot(2,3,1)
semilogy(NEPmin.KIDID,NEPmin.NEPGRquick,'sk','MarkerSize',4,'MarkerfaceColor','k');  hold on
xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
title('Phase readout')
subplot(2,3,3)
semilogy(NEPmin.KIDID,NEPmin.NEPGRquick,'sk','MarkerSize',4,'MarkerfaceColor','k');  hold on
xlabel('KID ID');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
title('R readout')
subplot(2,3,4)
semilogy(NEPmin.AluLength,NEPmin.NEPGRquick,'sk','MarkerSize',4,'MarkerfaceColor','k');  hold on
xlabel('Alu Length  (\mum)');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
title('Phase readout')
subplot(2,3,6)
semilogy(NEPmin.AluLength,NEPmin.NEPGRquick,'sk','MarkerSize',4,'MarkerfaceColor','k');  hold on
xlabel('Alu Length  (\mum)');ylabel('NEP (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
title('R readout')

for nm = 1:2:3
    subplot(2,3,nm)
    if opticalnepisthere == 1
        legend('optical NEP','Dark NEP','GR' )
    else
        legend('Dark NEP','GR' )
    end
end

subplot(2,3,2)
if opticalnepisthere == 1
    legend('dx/dP_{abs}','dx/dP_{dark}','10^{20} x dX/dN_{qp}','Location','best');
else
    legend('dx/dP_{dark}','10^{20} x dX/dN_{qp}','Location','best');
end
subplot(2,3,5)
if opticalnepisthere == 1
    legend('dx/dP_{abs}','dx/dP_{dark}','10^{20} x dX/dN_{qp}','Location','best');
else
    legend('dx/dP_{dark}','10^{20} x dX/dN_{qp}','Location','best');
end
MakeGoodFigure(15,9,12)



if opticalnepisthere == 1
    % optional debug figure only vs KID ID
    figure(2000)
    %optical data
    subplot(2,2,1)
    semilogy(OptNEPmin.KIDID,OptNEPmin.theta_Pabs,'or','MarkerSize',8,'MarkerfaceColor','r');  hold on
    
    subplot(2,2,2) %efficiency
    plot(OptNEPmin.KIDID,OptNEPmin.Efficiency,'ok','MarkerSize',8,'MarkerfaceColor','k');  hold on
    xlabel('KID ID');ylabel('optical Efficiency');grid on;
    
    subplot(2,2,3)%responsivity
    semilogy(OptNEPmin.KIDID,OptNEPmin.dthetadP,'or','MarkerSize',8);  hold on
    semilogy(OptNEPmin.KIDID,OptNEPmin.dthetadP./OptNEPmin.Efficiency,'or','MarkerSize',8,'MarkerfaceColor','r');  hold on
    subplot(2,2,4)%phase noise level
    semilogy(OptNEPmin.KIDID,OptNEPmin.phaserootNoiselevelfref,'or','MarkerSize',8,'MarkerfaceColor','r');  hold on
    
    % dark data
    subplot(2,2,1)
    semilogy(NEPmin.KIDID,NEPmin.theta,'ob','MarkerSize',6,'MarkerfaceColor','b');  hold on
    xlabel('KID ID');ylabel('NEP(\theta) (W /\surd Hz)');grid on;ylim([0.5e-20 1e-17]);
    legend('optical NEP(P_{abs})','Dark NEP')
    
    subplot(2,2,3)%responsivity NEPmin.Stheta_NEP
    semilogy(NEPmin.KIDID,NEPmin.dthetaPdark,'ob','MarkerSize',8,'MarkerfaceColor','b');  hold on
    xlabel('KID ID');ylabel('d\theta/dP');grid on;
    legend('d\theta dP','d\theta dP_{abs}','d\theta dP_{dark}')
    
    subplot(2,2,4)%Noise level
    semilogy(NEPmin.KIDID,NEPmin.Stheta_NEP,'ob','MarkerSize',8,'MarkerfaceColor','b');  hold on
    xlabel('KID ID');ylabel('\surd{S\theta} (1/Hz)');grid on;
    legend('Optical data','dark data')
    
MakeGoodFigure(9,9,12)
end
end