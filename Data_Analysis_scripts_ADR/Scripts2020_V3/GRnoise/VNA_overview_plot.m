figure
plot(VNAsweep.VarName1,VNAsweep.VarName2,'LineWidth',1)
fdesign = [4.6 4.7 4.8 5.1 5.2 5.3 5.7 5.8 5.9 6.2 6.3 6.4];
handleVisible = [{'on'},{'off'},{'off'}]
for i=1:3 
    xline(fdesign(i),'--','Color','red','LineWidth',1,'HandleVisibility',handleVisible{i})
end
for i=1:3 
    xline(fdesign(i+3),'--','Color','green','LineWidth',1,'HandleVisibility',handleVisible{i})
end
for i=1:3 
    xline(fdesign(i+6),'--','Color','blue','LineWidth',1,'HandleVisibility',handleVisible{i})
end
for i=1:3 
    xline(fdesign(i+9),'--','Color','magenta','LineWidth',1,'HandleVisibility',handleVisible{i})
end

xlim([4.5,6.5]);grid on;ylim([-40,5])
legend('VNA sweep','G1 design F_{0}','G2 design F_{0}','G3 design F_{0}','G4 design F_{0}')
title('VNA sweep')
xlabel('F [GHz]')
ylabel('|S_{21}| [GHz]')