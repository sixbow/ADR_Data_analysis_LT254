
function findrad_chip(nFG)
% same as main function, but siplified to be valid for the last (chip_
% stage

global FG  S21tot Fr S21tot_signal ddf 



FG(nFG).beamangle = atan(0.5*FG(1).D/FG(nFG).dist);     %beam angle entering: atan of TP def aperture and distance
FG(nFG).Bri = FG(nFG-1).Bri .* FG(nFG-1).S21;           %brilliance falling on this stage: Bri * S21 of previous filter
FG(nFG).S21 = ones(size(Fr));                       % Filter transmission = 0.5
FG(nFG).S21_signal = ones(size(Fr));                    % Filter transmission signal
S21tot = S21tot .* FG(nFG).S21;
S21tot_signal = S21tot_signal .* FG(nFG).S21_signal;

FG(nFG).Area = pi/4*(1e-3*FG(nFG).D)^2;                     %area current filter
FG(nFG).Etendue = pi*(1-cos(FG(nFG).beamangle)^2) * FG(nFG).Area; %Etendue - Omega * Area, Omega inc. reduction of filter area on large angle%
% note: solidangle2=2*pi*(1-cos(beamangle)); for uncorrected version
FG(nFG).Irad = FG(nFG).Bri * FG(nFG).Etendue;               %iradiance

% Power coupled from a previous Cardiff filter that gets hot: This makes FG(nFG).Irad_filt ~= 0 %
if ~isnan(FG(nFG).dist2F) %will calcualte only if there is a filter radiating. 
    if isfield(FG(nFG-1),'D') || isfield(FG(nFG),'Area')
        if isempty(FG(nFG).D) || isempty(FG(nFG).Area)
            error(['requested to calculate the power from a prev. hot filter, but for nFG: ' num2str(nFG-1) ' angle and etendue absent'])
        end
    end
    FG(nFG).filterangle = atan(0.5*FG(nFG-1).D/FG(nFG).dist2F);                     %angle to prev filter:
    FG(nFG).filterdilution = 1*pi*(1-cos(FG(nFG).filterangle)^2) / (1*pi);          %fractional TP from prev filter
    %The irradiance from emission prev filter = Irradiation absorbed (n-1) * etenude/2pi.%
    %The normalization with 2pi is needed because we have left the correct
    %normilzization of the brilliance. We assume al power is radiated in a half sphere., This is similar to SY filter dillution reasoning %
    FG(nFG).Irad_filt = FG(nFG-1).Iradabs * FG(nFG).filterdilution;                 %The irradiance from emission prev filter: Irad * solidangle/2pi
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

% calculating the power absorbed and irradiance absorbed 
FG(nFG).Iradabs = FG(nFG).Iradtot;  % iradiance absorbed
FG(nFG).Pabs = sum(FG(nFG).Iradabs) * ddf;                      % power absorbed in filterp

%This is not a cardiff filter, no filter heating
    FG(nFG).Tfilter = NaN;

end