clc 
clear all 
close all

dir_DCfiles = '../../Data_LT254_Sietse/LT254_DC_Chip13/';

Al_up_file = 'Al_up.csv';
Al_down_file = 'Al_down.csv';
%dir([dir_DCfiles Al_up_file])
DC_Al_up = csvread([dir_DCfiles Al_up_file])
DC_Al_down = csvread([dir_DCfiles Al_down_file])
hold on
yline(3.816,'--')
yline(3.816*0.9,'--')
yline(3.816*0.5,'--')
yline(3.816*0.1,'--')
plot(DC_Al_up(:,1),DC_Al_up(:,3),DC_Al_down(:,1),DC_Al_down(:,3))
hold off
xlim([1,1.4]);grid minor;ylim([0,4])



