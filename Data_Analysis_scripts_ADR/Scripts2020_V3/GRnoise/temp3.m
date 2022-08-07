clc
clear all 
close all 

a=[34.8 31.2 29 26.7 39.5];%dummy data
n=33;
[~,~,idx]=unique(round(abs(a-n)),'stable');
minVal=a(idx==1)