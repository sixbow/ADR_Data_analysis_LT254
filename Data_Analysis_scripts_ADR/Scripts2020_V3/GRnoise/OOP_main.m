%% Refactor OOP version!


%% (1) Single plots
o = Cfit(FFTsubsubdir,'NOISE_2D.mat')
[CTLS,gamma]= o.genfitsingle(1,5,10);
o = o.init_figax(2,2,'loglin');
o.plotsingle(2,2,1,5,10)









