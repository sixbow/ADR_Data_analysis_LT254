function [Tf,Pabs,Pemit_cold,Irademit_coldside,Iradabs] = getfilterT2(Iradin,Filter,Fr,Afilter,helpplot,SD,Shader)
% Iradin is the  irradiance, of the input radiation in W/Hz, brilliance * area * solidangle%
% Fr is the frequency array corresponding to the two previous vaklues, (Hz)%
% Filter is the ususal struct describibg a Cardiff filter
% Filter{nfi}(1,:) is transmission
% Filter{nfi}(2,:) is refelction
% Filter{nfi}(3,:) is absorption
% Shader is the same for  a cardiff shader, which has only row 1 (transmission)%
% Afilter is fiter area in m^2
% SD is the fractional solidangle between shader and filter (0..1)
%
% Tf = filter T (K)
% Irademit_filter emitted irradiance (W/Hz)
% Iradabs = absorebd irradiance
% Refls = 1D array of reflection ceofficient of the previous shader. 
% ReflsD is the fraction of reflected power comming back to the filter (typically solidangle/pi) 
%
% function calcualtes the temperature of a single sheet absorber shining
% into a 0K black body, loaded with SpecIr from one side
% valid for a flat surface, i..e TP = 2*pi*area, for radiation into two
% half spheres (not 4pi, because at 90 degree there is no emisison)


% constants
h = 6.6262e-34;		% J.s
c = 2.9998e8;		% m/s
k = 1.3806e-23;		% J/K
%fixed inputs
pol = 2;
ddf = Fr(2)-Fr(1);

% help plot
if nargin == 4
    helpplot = 0;
    SD = 0;Shader = [];
elseif nargin == 5
    SD = 0;Shader = [];
end

% anonymous fy giving the irradiance of the filter in a full sphere.
% the 2pi is because the filter is flat, not a sphere (then 4pi would be the
% total solidangle)
% Brilliance * emissivity 
Ira = @(T) (pol*h*Fr.^3)/(c^2).* 1./(exp(h*Fr./(k.*T))-1) .* Filter(3,:) * Afilter * 2 * pi;

%construct T vs P array
Ts = [2:2:310];
Pemit = zeros(length(Ts),1);
for n=1:length(Ts)
    Irademit = Ira(Ts(n));   % Irradiation of the warm filter (total from both sides)
    Pemit(n) = sum(Irademit) * ddf;
end

%power in: get filter T
Iradabs = Iradin .* Filter(3,:);  %Irradiation absorbed
Pabs = sum(Iradabs) * ddf ;
Tstart = interp1(Pemit,Ts,Pabs);        %T start temperature

%now iterative loop to determine Tf
% Input Irradiance to the filter: sum = input Irad
Iin = Iradin; clear Iradin  % Irad comming in through shader. This is the NET INPUT IRAD TO THE SYSTEM
Ir1r = zeros(size(Iin));     % Irad reflected from filter and shader comming back into filter
Ie1r = zeros(size(Iin));    % Irad emitted from filter and shader comming back into filter
% Output Irradiance from the filter
Ir1 = zeros(size(Iin));      % Irad reflected from filter from all power absorbed
Ie1 = zeros(size(Iin));     % Irad emitted from filter @ shader side
Ie2 = zeros(size(Iin));     % Irad emitted from filter @ cold side
% Irad leaving the system
Ir1o = zeros(size(Iin));     % Irad reflected from filter and transmitted out through the shader
Ie1o = zeros(size(Iin));    % Irad emitted from filter and transmitted out through the shader
Ie2o = zeros(size(Iin));    % Irad emitted from filter backside
Ir1a = zeros(size(Iin));     % Irad reflected from filter and absorbed before shadre
Ie1a = zeros(size(Iin));    % Irad emitted from filter and absorbed before shadre
Ie2 = zeros(size(Iin));     % Irad emitted from filter and transmitted out through the shader, already defined
Ito = zeros(size(Iin));      % Irad transmitted through the filter

Iradabs = (Iin + Ir1r + Ie1r) .* Filter(3,:);       % Iradiance absorbed in the filter
pf = 1;
if nargin == 7 %shader defined, full iterative calcul
    for ni = 1:100
        % power absorbed in the filter
        Ifin = (Iin + Ir1r + Ie1r); %net Irad on the filter. the p[refactor pf allows to reduce the steps
        Iradabs = Ifin .* Filter(3,:);          % Iradiance absorbed the filter
        Pabs = sum(Iradabs) * ddf ;
        % filter T
        Tf = interp1(Pemit,Ts,Pabs);                        %T filter where Pemit = Pabs
        
        % power leaving filter
        Ie1 = Ira(Tf)/2;
        Ie2 = Ira(Tf)/2;
        It2 = Ifin .* Filter(1,:);              % Iradiance transmitted through the filter
        Ir1 = Ifin .* Filter(2,:);              % Iradiance reflected from the filter
        
        %const check @ filter: Pin = Pout
        Pback_ori = sum(Ir1r + Ie1r) * ddf;
        Pfin_ori = sum(Ifin) * ddf ; %Ori input
        Pt2 = sum(It2) *ddf;
        Pr1 = sum(Ir1) *ddf;
        Pfout_ori = Pt2 + Pr1 + Pabs;
        %[Pfin_ori Pfout_ori Pr1 Pt2 Pabs]%ok
        
        % power comming back to the filter from Ir1 and Ie1
        Ir1r = Ir1 .* SD .* (1-Shader(1,:));   %filter and shader reflected
        Ie1r = Ie1 .* SD .* (1-Shader(1,:));  %filter emitted, shader reflected
        
        % power leaving the system 1: through the shader from Ir1 and Ie1
        Ie1o = Ie1 .* SD .* Shader(1,:);      % emission and reflection filter away through the shader
        Ir1o = Ir1 .* SD .* Shader(1,:);      % emission and reflection filter away through the shader
        % power leaving the system 2: absorbed between shader and filter%
        Ie1a = Ie1 .* (1-SD) ;      % emission and reflection filter away through the shader
        Ir1a = Ir1r .* (1-SD) ;      % emission and reflection filter away through the shader
        % power leaving the system 3: through the backside
        %Ie2o = Ie2; %just as a reminder
        %Ito = It2; %just as a reminder
        
        % checks @ filter level
        Pr1check = sum(Ir1r + Ir1o + Ir1a) * ddf; %reflected powers cons. check
        %[Pr1 Pr1check ]%ok
        Pe1 = sum(Ie1) * ddf; Pe1check = sum(Ie1r + Ie1o + Ie1a) * ddf; %emitted powers cons. check
        %[Pe1 Pe1check]%ok
        Pe = sum(Ie1 + Ie2 ) * ddf; % filter level Pabs and Pe check, should be close (not 1:1 dye to interp1)
        %[Pabs Pe]%about ok
        % power leaving system: 
        Pout = sum(Ie1o + Ir1o + Ie1a + Ir1a + Ie2 + It2) *ddf;
        % power into system
        Pin = sum(Iin) * ddf; %does not change, power from outside (for next loop)
        bla(ni,:) = [Pout Pin Tf];
        % new situation @ filter
        Pback = sum(Ir1r + Ie1r) * ddf; %backreflected power (for next loop)
        Pnextrun = Pback + Pin;
        %[Pfin_ori Pback Pnextrun] %ok
        
        if abs((Pout-Pin)/Pin) < 0.001 
            break
        else
            
        end
        
    end
else
    Tf = Tstart ;   %equilibrium when Pabs = Pemit
end



if helpplot == 1
    figure(1974)
   plot(1:ni,bla(1:ni,2));hold on
   plot(1:ni,bla(1:ni,1));
   plot(1:ni,(bla(1:ni,1) - bla(1:ni,2))./bla(1:ni,2));
   grid on
   xlabel('index');ylabel('power  (W)');axis tight
   legend('Pin','Pout')
end

%calcaul of the filter irradiance and so on uinder equilibrium
Irademit_coldside = Ira(Tf)/2;       %  Irradiation of the warm filter to the cold
Pemit_cold = sum(Irademit_coldside) * ddf ;
end