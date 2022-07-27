
function findrad(nFG,Filter,Tsource)
% main function of MAIN2_filtersAdvanced
% finds all power and filter T
% -- if  nFG == 1 system assumes this is the window, and that it is a BB
% with small reflection. So S21 == 1*(1-Refl) for power calcualtion on
% the stage (.S21), and given by S21(HDPE) for the signal 
% Calculates area of current filter and angle, TP to previous filter (Using .D and .dist)%
% -- if .CardiffFilter == 1, then it calculates the power absorbed in this
% filter, the irradiance absorbed and the temperature it will get, assuming a non-thermalized
% filter with a reflection of RC (0..1)
% -- if .dist2F is not empty, it will take the absorbed power of the previous
% filter = emitted power and calcualte the fractional solid angle using
% .dist2F to calcualte the TP and thereby tyhe p[ower radiated from the hot
% previous filter to the curren one.
%
% All length units are in mm!!!!
% % polyprop is a 2 row matrix with 1: transmission and 2 absoirption
global FG nPol S21tot Fr S21tot_signal ddf

h = 6.6262e-34;		% J.s
c = 2.9998e8;		% m/s
k = 1.3806e-23;		% J/K

%source brilliance for 2 polarizations!
pol = 2;
if nargin == 3
    Tsource = 293;
end
occupation=1./(exp(h*Fr/(k*Tsource))-1);
brilliance=(pol*h*Fr.^3)/(c^2).*occupation; %W/Hz/setrrad/m^2


if nFG == 1 %window
    FG(nFG).beamangle = pi/2;                               %beam angle entering in rad
    FG(nFG).refl = Filter{FG(nFG).Fi}(2,:);                 %Window reflection
    FG(nFG).Bri = brilliance;                               %brilliance falling on this stage (W/Hz/m^2/sterrad)
    FG(nFG).S21 = 1-FG(nFG).refl;                           %transmission for total power: window radiates of transmits 300K
    FG(nFG).S21_signal = Filter{FG(nFG).Fi}(1,:);           %Filter transmission for signal = window transmission
    S21tot = FG(nFG).S21;
    S21tot_signal = FG(nFG).S21_signal;
else %not window
    FG(nFG).beamangle = atan(0.5*FG(1).D/FG(nFG).dist);     %beam angle entering: atan of TP def aperture and distance
    FG(nFG).Bri = FG(nFG-1).Bri .* FG(nFG-1).S21;           %brilliance falling on this stage: Bri * S21 of previous filter
    if nFG ~= nPol %not polarizer
        FG(nFG).S21 = totalS21(Filter(FG(nFG).Fi));                  % Filter transmission total power
        FG(nFG).S21_signal = totalS21(Filter(FG(nFG).Fi));           % Filter transmission signal
        
    elseif nFG == nPol
        FG(nFG).S21 = 0.5*ones(size(Fr));                       % Filter transmission = 0.5
        FG(nFG).S21_signal = ones(size(Fr));                    % Filter transmission signal
    end
    S21tot = S21tot .* FG(nFG).S21;
    S21tot_signal = S21tot_signal .* FG(nFG).S21_signal;
end

FG(nFG).Area = pi/4*(1e-3*FG(nFG).D)^2;                     %area current filter in m^2
FG(nFG).Etendue = 1*pi*(1-cos(FG(nFG).beamangle)^2) * FG(nFG).Area; %Etendue - Omega * Area
FG(nFG).Irad = FG(nFG).Bri * FG(nFG).Etendue;               %iradiance W/Hz

% Power coupled from a previous Cardiff filter that gets hot: This makes FG(nFG).Irad_filt ~= 0 %
if ~isnan(FG(nFG).dist2F) %will calcualte only if there is a filter radiating. 
    if isfield(FG(nFG-1),'D') || isfield(FG(nFG),'Area')
        if isempty(FG(nFG).D) || isempty(FG(nFG).Area)
            error(['requested to calculate the power from a prev. hot filter, but for nFG: ' num2str(nFG-1) ' angle and etendue absent'])
        end
    end
    FG(nFG).filterangle = atan(0.5*FG(nFG-1).D/FG(nFG).dist2F);                     %angle to prev filter:
    FG(nFG).filterdilution = 1*pi*(1-cos(FG(nFG).filterangle)^2) / (1*pi);          %fractional TP from prev filter; 1pi%
    FG(nFG).Irad_filt = FG(nFG-1).Irademitcold * FG(nFG).filterdilution;            %The irradiance from emission prev filter: Irad * solidangle/2pi
else %no filter before
     FG(nFG).Irad_filt = 0;      %no radiation from hot filetr before this obne
     FG(nFG).filterangle = NaN;
     FG(nFG).filterdilution = NaN;
end

% getting powers
FG(nFG).Iradtot = FG(nFG).Irad_filt + FG(nFG).Irad;         % total irradiation
FG(nFG).Pdirect = sum(FG(nFG).Irad) * ddf ;                 % power entering this stage from 300K
FG(nFG).Pfilter = sum(FG(nFG).Irad_filt) * ddf ;            % power entering this stage from a hot filter before
FG(nFG).P = FG(nFG).Pdirect + FG(nFG).Pfilter;              % total power entering this stage

% calculating the power absorbed and irradiance absorbed in a Cardiff filter%
if FG(nFG).CardiffFilter == 1 %if this is a Cardiff filter
    hoek = atan(0.5*FG(nFG).D/FG(nFG).dist2shader);
    SD = 1*pi*(1-cos(hoek)^2) / (1*pi); 
    clear hoek;
    if ~isnan(FG(nFG).dist2shader) && length(FG(nFG).dist2shader) == 1 %single number: there is a shader before
        for ns=1:length(FG(nFG-1).Fi)
            if ns == 1
                Shadertot = Filter{FG(nFG-1).Fi(ns)};
            else
                Shadertot = Shadertot .* Filter{FG(nFG-1).Fi(ns)};
            end
        end
        [FG(nFG).Tfilter,FG(nFG).Pabs,FG(nFG).Pemitcold,FG(nFG).Irademitcold,FG(nFG).Iradabs] = ...
            getfilterT2(FG(nFG).Iradtot,Filter{FG(nFG).Fi},Fr,FG(nFG).Area,0,SD,Shadertot);
    else
        [FG(nFG).Tfilter,FG(nFG).Pabs,FG(nFG).Pemitcold,FG(nFG).Irademitcold,FG(nFG).Iradabs] = ...
            getfilterT2(FG(nFG).Iradtot,Filter{FG(nFG).Fi},Fr,FG(nFG).Area); %no shader before
    end 
else %This is not a cardiff filter, no filter heating
    FG(nFG).Tfilter = NaN;
    FG(nFG).Iradabs = zeros(size(FG(nFG).Iradtot));    %no radiation from hot filetr before this obne
    FG(nFG).Pabs = 0;
end

end
