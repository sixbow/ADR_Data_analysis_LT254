
function [TotalPbb,NEP,freq,filter]=BBradiator_mulitmode(Tbb, fstart, fstop ,transmission, Area, beamangle, delta, eta_pb, dF_GHz)
% 2019
% sub VI to get BB power
% calculates power emitted from blackbody at all Tbb temperatures
% between fstart and Fstop (so not outside the filter band)
% assumes 0 transmission outside filter band
%
% Input
% Tbb = T in K of the radiator. This can be an array.
% Fstart, Fstop are in THz
% transmission 0<=transmission<=1 is transmision of the optical chain
% OPTIONAL 
% delta is gap in J. Default 2.9818e-23J = 45 GHz =  Aluminium
% eta_pb is the pair breaking efficienyc, default = 0.4.. 
% dF_GHz            %integration steps in GHz. Defaultt = 1
% Area is receiving area in m^2
% beamangle = beam angle with respect to the normal in degrees (so total opening angle twice larger), in rad
% 
% Output
% TotalPbb = power in W in one mode, 1 polarization, transmitted
% NEP is struct witrh fields:
% NEP.g_r(p)=sqrt(sum(Flannigan*powerbla*delta/eta_pb));    %g-r NEP, only recombination
% NEP.poisson(p)=sqrt(sum(2*powerbla.*h.*freq));            %Poisson NEP
% NEP.wave(p)=sqrt(sum(wave));                              %NEP wave contribution
% NEP.totphoton=sqrt(NEP.g_r.^2+NEP.poisson.^2+NEP.wave.^2);%total KID NEP
% NEP.photonnoise=(NEP.g_r.^2+NEP.poisson.^2+NEP.wave.^2);  %total photon noise


%constants
Flannigan = 4; %4
debugplotje = 0; %1 plotrs a plotje
h = 6.6262e-34;		% J.s
c = 2.9998e8;		% m/s
k = 1.3806e-23;		% J/K
pol=1;                  %# of polarziation modes
Delta_GHz = 45.0;

if nargin == 6
    eta_pb=0.4;
    delta = Delta_GHz*1e9*h; % gap value in J
    df=1E9;             %integration steps
elseif nargin == 7
    eta_pb = 0.4;
    df=1E9;             %integration steps
elseif nargin == 8
    df=1E9;             %integration steps
elseif nargin == 9
    df=dF_GHz*1e9;             %integration steps
end

% making top hat F filter
fstart=fstart*1E12;%to Hz
fstop=fstop*1E12;
freq=fstart:df:fstop;

npts=length(freq);filter=zeros(1,npts);
if df*2 > fstop-fstart
    error('BBradiator_general: not emough resolution inside BBscript to calculate the requested BW')
end
for n=1:npts
    if freq(n)>=fstart && freq(n)<=fstop
        filter(n)=transmission;
    else
        filter(n)=0;
    end
end

% initialize arrays
N_Tbb=length(Tbb); %number of blackbody temperatures
NEP.g_r=zeros(1,N_Tbb); %GR noise
NEP.poisson=zeros(1,N_Tbb); %Poissonian photon noise
NEP.wave=zeros(1,N_Tbb); %Wavebunching photon noise
NEP.totphoton=zeros(1,N_Tbb); %Total noise due to photon fluctuations
TotalPbb=zeros(1,N_Tbb);
NEP.photonnoise = zeros(1,N_Tbb);

%solidangle
solidangle=1*pi*(1-cos(beamangle)^2); %including  effective rea reduction at large angles

for p=1:length(Tbb)
    %Photon Occupation Number (Bose-Einstein distribution):
    occupation=1./(exp(h*freq/(k*Tbb(p)))-1);
    
    %Planck brillance W/(m^2 sterrad Hz) [B_{\nu}(T_{bb})]:
    brilliance=(pol*h*freq.^3)/(c^2).*occupation;
    
    % total irradiation with fixed throughput but without filters in W/Hz
    % [B_{\nu}(T_{bb})*A\Omega]:
    % throughput,'lambda^2'
    Etendue=solidangle*Area+zeros(1,length(freq)); %array
    
    irradiation=brilliance.*Etendue;
    % total irradiation @ lens front in W/Hz
    % [B_{\nu}(T_{bb})*A\Omega*F_{\nu}]:
    filterirad=irradiation.*filter;%NB: This does not do anhthing as irradiation is also calcualetd only over filter BW
    
    % Power per frequency bin in W.
    % [B_{\nu}(T_{bb})*A\Omega*F_{\nu} d\nu]:
    powerbla=filterirad*df;
    
    % Total received power [W] Summation over power in all F bins
    % [\int_0^\infty B_{\nu}(T_{bb})*A\Omega*F_{\nu} d\nu]:
    TotalPbb(p)=sum(powerbla);
    
    % Calculation of the various NEP's
    NEP.g_r(p)=sqrt(sum(Flannigan*powerbla*delta/eta_pb)); %g-r noise, only recombination
    NEP.poisson(p)=sqrt(sum(2*powerbla.*h.*freq)); %Poisson noise term
    
    lambda2=(c./freq).^2; %lamdba^2
    % wave bunching: 2Phf        *             opt eff          *    B
    wave=(2*powerbla.*h.*freq).*(filter.*Etendue./lambda2).*occupation;
    NEP.wave(p)=sqrt(sum(wave)); %wave contribution
    NEP.filterirad(p,:) = filterirad; %This is identicalk to irradiation because no data ouside filter F range
    NEP.F = freq;
end

if debugplotje ==1
    figure(1)
    subplot(1,2,1)
    semilogy(Tbb,TotalPbb,'.-r','MarkerSize',8);hold on
    ylabel('Power in W');legend('filtered')
    xlabel('T (K)')
    grid on
    
    subplot(1,2,2)
    plot(freq,filter)
    ylabel('Transmission')
    xlabel('Frequency (Hz)')
    ylim([0 1.1])
    axis tight;
    xlim([0 1e12])
    grid on
end

NEP.totphoton=sqrt(NEP.g_r.^2+NEP.poisson.^2+NEP.wave.^2); %total KID NEP
NEP.photonnoise=NEP.g_r.^2+NEP.poisson.^2+NEP.wave.^2; %total photon noise
end

