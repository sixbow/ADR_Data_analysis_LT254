function [Tf,Irademit_filter,Iradabs] = getfilterT(Iradin,eps,Fr,Afilter,Refl,helpplot)
% Iradin is the  irradiance, of the input radiation in W/Hz, brilliance * area * solidangle%
% eps is the spectral absorption (or emissivity) (1/Hz)
% Fr is the frequency array corresponding to the two previous vaklues, (Hz)%
% Tstart = optional is start filter T (K)
% Afilter is the filter area in m^2
% Refl is the reflection coefficient of the filter
%
% Tf = filter T (K)
% Irademit_filter emitted irradiance (W/Hz)
% Iradabs = absorebd irradiance
%
% function calcualtes the temperature of a single sheet absorber shining
% into a 0K black body, loaded with SpecIr from one side
% valid for a flat surface, i..e TP = 2*pi*area, for radiation into two
% half spheres (not 4pi, sive at 90 degree there is noe misison)


% constants
h = 6.6262e-34;		% J.s
c = 2.9998e8;		% m/s
k = 1.3806e-23;		% J/K
%fixed inputs
pol = 2;
ddf = Fr(2)-Fr(1);

% help plot
if nargin == 5
    helpplot = 0;
end

%power in
Iradabs = (1 - Refl) * Iradin .* eps;  %correcting for filter absorption spectrum
Pabs = sum(Iradabs) * ddf ;

%construct T vs P array
Ts = [2:2:50 55:5:310];
Pemit = zeros(length(Ts),1);
for n=1:length(Ts)
    occupation=1./(exp(h*Fr./(k.*Ts(n)))-1);
    Bri=(pol*h*Fr.^3)/(c^2).*occupation;            %brilliance W/Hz/m^2/solidangle
    Irademit = (1 - Refl) * Bri .* eps .* Afilter * 2* pi;   % Irradiation of the warm filter (total from both sides)
    Pemit(n) = sum(Irademit) * ddf;
end

Tf = interp1(Pemit,Ts,Pabs) ;   %equilibrium when Pabs = Pemit

if helpplot == 1
   plot(Ts,Pemit);hold on
   plot(Tf, Pabs,'or','markerfacecolor','r');grid on
   xlabel('T  (K)');ylabel('Power  (W)');axis tight
end

%calcaul of the filter irradiance and so on uinder equilibrium
occupation=1./(exp(h*Fr./(k.*Tf))-1);
Bri=(pol*h*Fr.^3)/(c^2).*occupation;            %brilliance W/Hz/m^2/solidangle
Irademit_filter = (1 - Refl) *  Bri .* eps .* Afilter * 2 * pi;       % total Irradiation of the warm filter
Pemitcheck = sum(Irademit_filter) * ddf ;
end