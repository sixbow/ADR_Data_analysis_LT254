close all;clear variables
% see als
%addpath('subscripts');

addpath(genpath('SubfilesGefilter'));
c = 2.9998e8;		% m/s
%set F range


[Filter,name,Fr] = Getfilters;
R = 2e-2;%window radius
Area = pi*R^2;
%name is a cell aray with short filter disriptions
%Filer is the cell array with filter data

%% design the filter stack
Fi = [18,17,19,6,24,23]; %this defines he filters used

%%  % FINAL FIGURES

set(0,'DefaultLegendAutoUpdate','off')
Fend = 1; %in THz

figure(1)
% 2 LPF @ holder: Baseline
T = 'resultsFinal/DeshimaFinal2LPFHolder';

Frange_inband = find(Fr > 230e9 & Fr< 250e9);   %in band power range (indices in Fr array)
figure(1)
[S21,Power] = PlotSystem(Fr,Area,Filter,name,Fi,Frange_inband);
subplot(3,1,1);xlim([0.1 Fend])
subplot(3,1,2);xlim([0.1 Fend])
plot(Fr/1e9,S21,'g','Linewidth',2);
subplot(3,1,3);ylim([1e-25 1e-15]);xlim([0.1 Fend]*1e12);
MakeGoodFigureBB(15,20,12,'Deshima_ADR_240');

%% data save
%save('FilterT.mat','Fr','S21')
outputdata = [Fr ; S21]';
dlmwrite(['Deshima_ADR_240_F_S21.dat'], outputdata,'delimiter','\t','precision',6);
outputdata2 = [0.01*Fr/c ; S21]'; %from F(Hz) to cm^-1
dlmwrite(['Deshima_ADR_240.dat'], outputdata2,'delimiter','\t','precision',6);
rmpath(genpath('Subfiles'));

