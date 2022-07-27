
h = 6.6262e-34;		% J.s
c = 2.9998e8;		% m/s
k = 1.3806e-23;		% J/K
%Define frequency range (for everything)
Fr = 0.01:0.01:150; %@ 150 THz (5e3cm-1) 300K radiation is negligible 

%source T
Tsource = 293;

Fr = Fr * 1e12;
ddf = Fr(2) - Fr(1);

% get filter specs
[Filter,name] = Getfilters(Fr);
ReradFreq = 400*100 *c;% at about 400 cm-1 filters will become emissive; and I assume they will also leak a bit.
Reradi = Fr > ReradFreq;

%source brilliance for 2 polarizations!
pol = 2;
occupation=1./(exp(h*Fr/(k*Tsource))-1);
brilliance=(pol*h*Fr.^3)/(c^2).*occupation;

% define frequency rage where filters become black bodies
% here I will assume that the Cardiff filters radiate all rpower F>ReradFreq in 2pi sterrad. Radiation at lower powers is treated as normal 
%ReradFreq = 400*100 *c;% at about 400 cm-1 filters will become emissive; and I assume they will also leak a bit.
%Reradi = Fr > ReradFreq;
%These values are known fron getfilter

% define filter groups
% window
nFG = 1; %
FG(nFG).Fi = 1;                         %index of 'Filter' that defines the filter properties
FG(nFG).S21 = Filter{FG(nFG).Fi}(1,:);  %Filter transmission
FG(nFG).refl = Filter{FG(nFG).Fi}(2,:); %Window reflection
FG(nFG).D = 80;                         %mm. diameter
FG(nFG).beamangle = 90;                 %beam angle entering
FG(nFG).Bri = brilliance;               %brilliance falling on this stage
FG(nFG).Area = pi/4*(1e-3*FG(nFG).D)^2;      %area current filter
FG(nFG).Etendue = 1*pi*(1-cos(FG(nFG).beamangle)^2) * FG(nFG).Area;     %Etendue - Omega * Area
FG(nFG).Irad = FG(nFG).Bri * FG(nFG).Etendue;   %iradiance
FG(nFG).P = sum(FG(nFG).Irad) * ddf;    %power entering this stage
FG(nFG).Pabs = 0;    %power absorbed

% Shaders at 300, 77 and 4K. Defined with entrance aperture
nFG = 2; %
FG(nFG).D = 70; %mm. diameter
FG(nFG).Fi = [2,  2,    2,    3,   3];                      %index of 'Filter' that defines the filter properties
FG(nFG).S21 = totalS21(Filter(FG(nFG).Fi));                 %Filter transmission
FG(nFG).dist = 340-327.6;                                   %mm. distance to TP defining aperture (=window)
FG(nFG).beamangle = atan(0.5*FG(1).D/FG(nFG).dist);      %beam angle entering: atan of TP def aperture and distance
FG(nFG).Area = pi/4*(1e-3*FG(nFG).D)^2;                          %area current filter
FG(nFG).Bri = brilliance .* (1-FG(nFG-1).refl);                %brilliance falling on this stage: 1-Refl(HDPE)
FG(nFG).Etendue = 1*pi*(1-cos(FG(nFG).beamangle)^2) * FG(nFG).Area;     %Etendue - Omega * Area
FG(nFG).Irad = FG(nFG).Bri * FG(nFG).Etendue;   %iradiance
FG(nFG).P = sum(FG(nFG).Irad) * ddf;    %power entering this stage
FG(nFG).Pabs = 0;    %power absorbed

% 4K LPF (first one). The angle to this filter is given by the window and
% filter distance
nFG = 3; %
FG(nFG).D = 37; %mm. diameter
FG(nFG).Fi = [4];                                           %index of 'Filter' that defines the filter properties
FG(nFG).S21 = totalS21(Filter(FG(nFG).Fi));                 %Filter transmission
FG(nFG).dist = 340-221.5;                                   %mm. distance to TP defining aperture (=window)
FG(nFG).beamangle = atan(0.5*FG(1).D/FG(nFG).dist);      %beam angle entering
FG(nFG).Area = pi/4*(1e-3*FG(nFG).D)^2;                          %area current filter
FG(nFG).Bri = FG(nFG-1).Bri .* FG(nFG-1).S21;               %brilliance falling on this stage: Bri prev stage * S21 prev stage
FG(nFG).Etendue = 1*pi*(1-cos(FG(nFG).beamangle)^2) * FG(nFG).Area;     %Etendue - Omega * Area
FG(nFG).Irad = FG(nFG).Bri * FG(nFG).Etendue;   %iradiance
FG(nFG).P = sum(FG(nFG).Irad) * ddf;    %power entering this stage
FG(nFG).Pabs = sum(FG(nFG).Irad(Reradi)) * ddf;    %power absorbed

% 4K LPF (2nd one). The angle to this filter is given by the window and
% filter distance
nFG = 4; %
FG(nFG).D = 37; %mm. diameter
FG(nFG).Fi = [5];                                           %index of 'Filter' that defines the filter properties
FG(nFG).S21 = totalS21(Filter(FG(nFG).Fi));                 %Filter transmission
FG(nFG).dist = 340-185.9;                                   %mm. distance to TP defining aperture (=window)
FG(nFG).beamangle = atan(0.5*FG(1).D/FG(nFG).dist);      %beam angle entering
FG(nFG).Area = pi/4*(1e-3*FG(nFG).D)^2;                          %area current filter
FG(nFG).Bri = FG(nFG-1).Bri .* FG(nFG-1).S21;               %brilliance falling on this stage: Bri prev stage * S21 prev stage
FG(nFG).Etendue = 1*pi*(1-cos(FG(nFG).beamangle)^2) * FG(nFG).Area;     %Etendue - Omega * Area
FG(nFG).Irad = FG(nFG).Bri * FG(nFG).Etendue;   %iradiance
FG(nFG).P = sum(FG(nFG).Irad) * ddf;    %power entering this stage
FG(nFG).Pabs = sum(FG(nFG).Irad(Reradi)) * ddf;    %power absorbed

% %% 4K LPF (pupil). TP defined by  4K filter  @optics box entrance, and defined in reverse. this TP is conserved for the rest of the system%
nFG = 5; %
FG(nFG).D = 47; %mm. diameter
FG(nFG).Fi = [6];                                           %index of 'Filter' that defines the filter properties
FG(nFG).S21 = totalS21(Filter(FG(nFG).Fi));                 %Filter transmission
FG(nFG).dist = 185.9-54.9;%221.5-54.9;                                   %mm. distance to TP defining aperture (=lyot stop)
FG(nFG).beamangle = atan(0.5*FG(3).D/FG(nFG).dist);      %beam angle entering
FG(nFG).Area = pi/4*(1e-3*FG(nFG).D)^2;                          %area current filter 
FG(nFG).Bri = FG(nFG-1).Bri .* FG(nFG-1).S21;               %brilliance falling on this stage: Bri prev stage * S21 prev stage
FG(nFG).Etendue = 1*pi*(1-cos(FG(nFG).beamangle)^2) * FG(nFG-1).Area;     %Etendue - Omega * Area, NB: area prev filter taken
FG(nFG).Irad = FG(nFG).Bri * FG(nFG).Etendue;   %iradiance for low frequencies\
FG(nFG).P = sum(FG(nFG).Irad) * ddf;    %power entering this stage
FG(nFG).Pabs = sum(FG(nFG).Irad(Reradi)) * ddf;    %power absorbed

% Polarizer@ end 4K box. The angle to this filter is given by 4K filter and total distance%
nFG = 6; %
FG(nFG).beamangle = NaN;      %beam angle entering
FG(nFG).D = 15; %mm. diameter
FG(nFG).S21 = 0.5*ones(size(Filter{1}(1,:)));                    %Filter transmission, =1/2
FG(nFG).Bri = FG(nFG-1).Bri .* FG(nFG-1).S21;               %brilliance falling on this stage: Bri prev stage * S21 prev stage
FG(nFG).Etendue = FG(5).Etendue ;
FG(nFG).Irad = FG(nFG).Bri * FG(nFG).Etendue;   %iradiance
FG(nFG).P = sum(FG(nFG).Irad) * ddf;    %power entering this stage
FG(nFG).Pabs = sum(FG(nFG).Irad(Reradi)) * ddf;    %power absorbed

% 100 mK holder. The angle to this filter is given by 4K filter and total distance%
nFG = 7; %
FG(nFG).beamangle = NaN;      %beam angle entering
FG(nFG).D = 12; %mm. diameter
FG(nFG).Fi = [7,8];                                           %index of 'Filter' that defines the filter properties
FG(nFG).S21 = totalS21(Filter(FG(nFG).Fi));                 %Filter transmission
FG(nFG).Bri = FG(nFG-1).Bri .* FG(nFG-1).S21;               %brilliance falling on this stage: Bri prev stage * S21 prev stage
FG(nFG).Etendue = FG(5).Etendue ;
FG(nFG).Irad = FG(nFG).Bri * FG(nFG).Etendue;   %iradiance
FG(nFG).P = sum(FG(nFG).Irad) * ddf;    %power entering this stage
FG(nFG).Pabs = sum(FG(nFG).Irad(Reradi)) * ddf;    %power absorbed
% Chip
nFG = 8; %
FG(nFG).beamangle = NaN;      %beam angle entering
FG(nFG).D = 10; %mm. diameter
FG(nFG).Bri = FG(nFG-1).Bri .* FG(nFG-1).S21;               %brilliance falling on this stage: Bri prev stage * S21 prev stage
FG(nFG).Etendue = FG(5).Etendue ;
FG(nFG).Irad = FG(nFG).Bri * FG(nFG).Etendue;   %iradiance
FG(nFG).P = sum(FG(nFG).Irad) * ddf;    %power entering this stage
FG(nFG).Pabs = sum(FG(nFG).Irad(Reradi)) * ddf;    %power absorbed

%%
kleur = colormapJetJB(nFG);
close all
figure(1)
subplot(1,2,1)
for n=1:nFG
    loglog(Fr,FG(n).Bri,'color',kleur(n,:));hold on
end
ylabel('Brilliance (W m^{-2} sterrad^{-1} Hz^{-1}))');grid on;
ylim([1e-25 1e-11])
subplot(1,2,2)
for n=1:nFG
    loglog(Fr,FG(n).Irad,'color',kleur(n,:));hold on
end
ylabel('Filtered Irradiance (W/Hz)')
grid on;ylim([1e-25 1e-11])
MakeGoodFigure(14,7,12)

figure(2)
subplot(2,2,1)
for n=1:nFG
    plot(n,FG(n).D,'o','color',kleur(n,:),'markerfacecolor',kleur(n,:));hold on
    ylabel('Diameter  (mm)')
end
grid on;
subplot(2,2,2)
for n=1:nFG
    semilogy(n,FG(n).beamangle,'o','color',kleur(n,:),'markerfacecolor',kleur(n,:));hold on
    ylabel('beamangle  (rad)')
end
grid on;
subplot(2,2,3)
for n=1:nFG
    semilogy(n,FG(n).Etendue,'o','color',kleur(n,:),'markerfacecolor',kleur(n,:));hold on
    ylabel('Etendue  (m^2 sterrad)')
end
grid on;
legend(['window','5 shaders 50K-4K',name(4:6),'Polarizer',[name{7} '+' name{8}], 'chip'])

subplot(2,2,4)
for n=1:nFG
    semilogy(n,FG(n).P,'o','color',kleur(n,:),'markerfacecolor',kleur(n,:));hold on
    semilogy(n,FG(n).Pabs,'o','color',kleur(n,:));
    ylabel('P({entering} (solid) and P_{absorbed} (open)  (W)');
end
ylim([1e-10 2])
grid on;
MakeGoodFigure(14,10,12)