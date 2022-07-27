function [Filter,name,Fr] = Getfilters(shaderleak,RC,z)
% Fr is frequency array in Hz
% Filter is a cell array with all filter transmissions etc.
% can run stand alone
% all Cardiff filters are patched at low F end with the lowest measured
% vakue. At the high end, they are set to Highfleak as long as they do not have
% data
%Cardiff shadres are patched with the lowets and highest datanumbers
%filter heating is defined using polypropylene absorption spectrum
% Filter{nF} is array of struct with fiulter infor
% namd is array of strings with names
% opolyprop is 2 rows, S21 and absorption of polypropylene, F range given
% by Fr
% Ts is Tempeature and P emitted by polypropylene per m and 2piu sterrad
%
% All filterdata is in ../TestedMatlabScripts/BlackBodies/filterfiles/'
%
% INPUTS all optional
% shaderleak is the transmission of the sahder at high F, default = 0 
% RC 0..1 is the refecrtion of a mesh filter at high F, default = 0
% z is the polypropylene thickness of the filters iof QMC. Defaul = 0.5mm
%
% polyprop is a 2 row matrix with 1: transmission and 2 absoirption
% PolyBri_Ts(1,:) is power/area/solidangle absorbed from a source of
% temperature Ts in polypropylene
% PolyBri_Ts(2,:) is associated source T
%
% Filter{nfi}(1,:) is transmission
% Filter{nfi}(2,:) is refelction
% Filter{nfi}(3,:) is absorption

FL = '/Users/jochem/ownCloud/KID/Technology_Experiment/TestedMatlabScripts/BlackBodies/filterfiles/';
if nargin == 0          %default settings for leaking shaders
    shaderleak = 0.2;   %0.2
    RC = 0.0;           %0
    z = 0.5e-3;         %0.5 mm thick filters (absorption calcul)
elseif nargin == 2
    z = 0.5e-3; %0.5 mm thick filters (absorption calcul)
end

helpplot = 0;

Fr = 0.08:0.002:12; %
%Fr = 0.1:0.3:10; %@ 150 THz (5e3cm-1) 300K radiation is negligible
Fr = Fr * 1e12;


h = 6.6262e-34;		% J.s
c = 2.9998e8;		% m/s
k = 1.3806e-23;		% J/K

minw = Fr(1)/c;
ReradFreq = 400*100 *c;% at about 400 cm-1 (12 THz) filters will become emissive; and I assume they will also leak a bit.
Reradi = Fr > ReradFreq;
Highfleak = 0.01;%high F leakage (transmission) to make things not too good
% Note that the coefficient RC defining excess filter reflection will only apply for this F range as well!

%read absorption data valid for all Cardiff filters. Taken from Tucker et
%al. col 1= wavenumber cm-1, col 2 absorption in Neper, checked to be Np/cm. ?
% defL: ?A/Ao = exp(-alpha_Np z)
polyprop_abs = csvread([FL 'absorption_filters.csv']); %wwavenumber, nepers/cm (i.e absortpion coeff alpha, )
polyprop_abs(:,3) = exp(-2*polyprop_abs(:,2)*z*100); %transmission due to PP absorption, z converted from m to cm, factor 2 from V to P
polyprop(1,:) = interp1( polyprop_abs(:,1)*100*c ,polyprop_abs(:,3),Fr); %polypropanol transmission due to absorption for Fr
polyprop(1,isnan(polyprop(1,:))) = 1; %assume transmisison is 1 when no data
polyprop(2,:) = 1 - polyprop(1,:); % power  absorption in the polyproanol vs Fr

% Filter 1: 8 mm thick HDPE window, from Stephen's king cryostat models
% measured data + interpolation based upon n=1.52 and dat from Lamb International Journal of lnfi'ared and Millimeter Waves, Vol, 17, No. 19., 1996
nfi = 1;
thickness   = 0.008;%m
reflens     = ((1-1.52)/(1+1.52))^2;
tandelta    = 4.805e-4; %% tan delta, measured Biorat
tan2delta   = 1e-8;%2.893e-8; %% tan delta, measured Biorat. I use 1e-8 as this fits the tail of the data better
load([FL 'filtersKing.mat'],'HDPEpup_clean')
AbsHDPE     = exp(-thickness*2*pi*1.52*(tandelta*Fr/c + tan2delta*(Fr/c).^2)).*(1-2*reflens);       %  transmisison model (incl refl.)
M_AbsHDPE   = interp1(HDPEpup_clean(:,1)*100*c,HDPEpup_clean(:,2)/100,Fr);                          % Measured transmission
M_AbsHDPE(M_AbsHDPE >1)=1;  %catching T > 1
HDPE        = M_AbsHDPE;
HDPE(isnan(HDPE)) = AbsHDPE(isnan(HDPE));
HDPE(HDPE<1e-3) = 1e-3;
Filter{nfi}(1,:) = HDPE;                              %Transmission (incl reflections): S21 = (1-abs)*(1-2*reflens)
Filter{nfi}(2,:) = ones(size(Fr))*reflens;            %reflection
Filter{nfi}(3,:) = 1 - HDPE/(1-2*reflens);            %absorption (or emission): 1-absorption = S21/(1-2*Reflens)
Filter{nfi}(4,:) = HDPE/(1-2*reflens);                %perfect HDPE transmission
name{nfi} = '8 mm HDPE';
clear tandelta tan2delta AbsHDPE M_AbsHDPE HDPE reflens thickness HDPEpup_clean
nfi =nfi + 1;

% Get Deshima1 filters
% Shaders
D = [FL 'Deshima/'];
%filter 2: C15
load('filtersKing.mat','C15')
Filter{nfi}(1,:) = interp1([Fr(1) C15(:,1)'*100*c Fr(end)], [1 C15(:,2)' C15(end,2)],Fr);
Filter{nfi}(1,Filter{nfi}(1,:) > 1) = 1;
Filter{nfi}(1,Filter{nfi}(1,:) < shaderleak) = shaderleak;
name{nfi} = 'C15';clear C15;
nfi = nfi + 1;
%filter 3: C30
load('filtersKing.mat','C30')
Filter{nfi}(1,:) = interp1([Fr(1) C30(:,1)'*100*c Fr(end)], [1 C30(:,2)' C30(end,2)],Fr);
Filter{nfi}(1,Filter{nfi}(1,:) > 1) = 1;
if helpplot == 1
    figure(12)
    loglog(Fr, Filter{nfi}(1,:));hold on;grid on
end
Filter{nfi}(1,Filter{nfi}(1,:) < shaderleak) = shaderleak;
name{nfi} = 'C30';clear C30;
if helpplot == 1
    loglog(Fr, Filter{nfi}(1,:));
    legend('Raw','leak adjusted ');xlabel('Frequency   (Hz)');ylabel('S21')
    MakeGoodFigure(15,6,14,'filtercheck/Shaderexample')
end
nfi = nfi + 1;

%filter 4: 2 THz LPF 
[F,S21] = ReadaFilter([D 'K2280_66cm-1_LPF.dat'],1);name{nfi} = 'K2280 2 THz LPF';%40K
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

%filter 5: 1.5 THz LPF
[F,S21] = ReadaFilter([D 'W2206_58cm-1_LPF.dat'],1);name{nfi} = 'W2206 1.75 THz LPF';%4K
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

%filter 6: 1.1 THz LPF
[F,S21] = ReadaFilter([D 'K2325_36cm-1_LPF.dat'],1);name{nfi} = 'K2325 1.1 THz LPF';%4K pupil
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

%% filter 7: 1.0 THz LPF (THIS IS THE EXAMPLE TO PLOT)
[F,S21] = ReadaFilter([D 'K1785_33cm-1_LPF.dat'],1);name{nfi} = 'K1785 1 THz LPF';%100 mK 1Thz LPF
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC, helpplot);
nfi = nfi + 1;
%%   
%filter 8: 350 GHz BPF
[F,S21] = ReadaFilter([D 'K1817_BPFDeshima.dat'],1);name{nfi} = 'K1817 350 GHz BPF';%100 mK 350 LPF
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

%get beast filters
B = [FL 'Beast/'];
%filter 9: RTMLI (1 layer)
load([B 'RTMLI.mat']);name{nfi} = 'RTMLI';%40K
F = RTMLI(:,1)';S21 = RTMLI(:,2)';clear RTMLI
Filter{nfi}(1,:) = abs(interp1([Fr(1) F ], [S21(1) S21], Fr));
Filter{nfi}(1,isnan(Filter{nfi}(1,:)) & Reradi) = Highfleak;Filter{nfi}(1,isnan(Filter{nfi}(1,:)) & ~Reradi) = abs(S21(end));
clear F S21;nfi = nfi + 1;

%filter 10: Gore Tex foam
load('filtersKing.mat','goretexfoam')%data in cm-1: *100*c to get Hz
Filter{nfi}(1,:) = interp1([Fr(1) goretexfoam(:,1)'*100*c Fr(end)], [1 goretexfoam(:,2)'/100 goretexfoam(end,2)/100],Fr);%NB: steve's data is in reverse F
Filter{nfi}(1,Filter{nfi}(1,:) > 1) = 1;Filter{nfi}(1,:) = abs(Filter{nfi}(1,:));
name{nfi} = 'goretexfoam';clear goretexfoam;
nfi = nfi + 1;

%filter 11: LPF
[F,S21] = ReadaFilter([B 'K2329.txt'],1);name{nfi} = 'K2329 1.5 THz LPF';
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

%Filter 12 LPF
[F,S21] = ReadaFilter([B 'W963.txt'],1);name{nfi} = 'W963 1 THz LPF';
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

%filter 13 7 RTMLI layers
Filter{nfi}(1,:) = Filter{9}(1,:).^7;name{nfi} = '7xRTMLI';%40K
Filter{nfi}(1,Filter{nfi}(1,:) < 0.005) = 0.5*0.01;%limiting how good these are
nfi = nfi + 1;
%filter 14 3 RTMLI layers
Filter{nfi}(1,:) = Filter{9}(1,:).^3;name{nfi} = '3xRTMLI';%40K
Filter{nfi}(1,Filter{nfi}(1,:) < 0.01) = 0.01;%limiting how good these are
nfi = nfi + 1;
%filter 15 perfect HDPE
Filter{nfi}(1,:) = Filter{1}(4,:);name{nfi} = 'AR coated HDPE';
Filter{nfi}(2,:) = ones(size(Fr))*0;            %reflection
nfi = nfi + 1;
%filter 16 Detector lens
Filter{nfi}(1,:) = ones(size(Filter{nfi-1}(1,:)));name{nfi} = 'detector';
nfi = nfi + 1;

%get ADR filters
AD = [FL];
%filter 17 B386
[F,S21] = ReadaFilter([AD 'B386 18cm LPE.dat'],1);name{nfi} = 'ADR B386 500 GHz LPF';% ADR 350 14cm-1 LPE  located in /testedmatlabscript 
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

%filter 18 W969
[F,S21] = ReadaFilter([AD 'W969 37cmLPESCUBAII.txt'],1);name{nfi} = 'ADR W969 1.2 THz LPF BB';%located in /testedmatlabscript 
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

%filter 19 W1052
[F,S21] = ReadaFilter([AD 'W1052 14cm LPE SCUBAII.dat'],1);name{nfi} = 'ADR W1052 400 GHz LPF' ;% ADR 350 14cm-1 LPE located in /testedmatlabscript 
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

%filter 20 W1275
[F,S21] = ReadaFilter([AD 'W1275_350GHz.dat'],1);name{nfi} = 'W1275 350 GHz BPF';%located in /testedmatlabscript 

Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

% Final Deshima2 filters
%Filter 21 New Deshima 630 GHz LPF
xldat = xlsread([D 'Deshima2Filters.xlsx'],'K2844','A1:B491');name{nfi} = 'K2844 21cm-1 with ARC';%located in /testedmatlabscript 
F=xldat(:,1)';  %
S21=xldat(:,2)';
F = F*100*c;
clear xldat;
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

%Filter 22 New Deshima 550 GHz LPF
xldat = xlsread([D 'Deshima2Filters.xlsx'],'K2845','A1:B491');name{nfi} = 'K2845 18cm-1 with ARC';%located in /testedmatlabscript 
F=xldat(:,1)';  %
S21=xldat(:,2)';
F = F*100*c;
clear xldat;
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

%Filter 23 New Deshima 450 LPF 15cm-1 LPE May 2016.xls.
xldat = xlsread([D 'Deshima2Filters.xlsx'],'K2846','A1:B589');name{nfi} = 'K2846 15cm-1 LPE';%located in /testedmatlabscript 
F=xldat(:,1)';  %
S21=xldat(:,2)';
F = F*100*c;
clear xldat;
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

%Filter 24 New Deshima 240 GHz BPF.
xldat = xlsread([D 'Deshima2Filters.xlsx'],'K2434','A1:B197');name{nfi} = 'K2434 240GHz BPF';%located in /testedmatlabscript 
F=xldat(:,1)';  %
S21=xldat(:,2)';
F = F*100*c;
clear xldat;
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

%Filter 25 New Deshima 350 GHz BPF.
xldat = xlsread([D 'Deshima2Filters.xlsx'],'K2584','A1:B354');name{nfi} = 'K2584 350GHz BPF';%located in /testedmatlabscript 
F=xldat(:,1)';  %
S21=xldat(:,2)';
F = F*100*c;
clear xldat;
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;

%Filter 26 New Deshima 650 GHz BPF.
xldat = xlsread([D 'Deshima2Filters.xlsx'],'K2136','A1:B589');name{nfi} = 'K2136 650GHz BPF';%located in /testedmatlabscript 
F=xldat(:,1)';  %
S21=xldat(:,2)';
F = F*100*c;
clear xldat;
Filter{nfi} = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC);
nfi = nfi + 1;


end
function Filter = maakfilter(Fr,F,S21, Reradi, polyprop, Highfleak, RC, helpplot)
if nargin == 7
    helpplot = 0;
end
%1 transmission
%2 reflection
%3 absorption

%create S21 (1) from data provided (and patch)%
Filter(1,:) = abs(interp1([Fr(1) F ], [S21(1) S21], Fr));  %patching to lowest freq + interpolate
RTT2 = isnan(Filter(1,:)) & ~Reradi;               %medium frequency range (F < 12 THz) has no data => patch last datapoint
RTT = isnan(Filter(1,:)) & Reradi;                 %high F > F.12 THz 
Filter(1,RTT2) = abs(S21(end));
Filter(1,RTT) = Highfleak;
%absorption (3)
Filter(3,:) = polyprop(2,:);                       % absortpion where there is data
tolarge = Filter(1,:) + Filter(2,:) > 1;      % catching where S21 + abs > 1
Filter(3,tolarge) = 1 - Filter(1,tolarge);    % abs = S21 - 1
% define reflection (2): S21 = 1-S11 -Abs 
Filter(2,:) = 1 - Filter(1,:)-Filter(3,:);
higherR = Filter(2,:) < RC & RTT;                 % indices where the reflection < RC AND at F > 12 THz
Filter(2,higherR) = RC;
% re-define absorption (3): we have made the reflectrion large wrt porev def of
% absorption: S21 = 1 - S11 - Abs so Abs = 1-S11-S21
Filter(3,higherR) = 1 - Filter(2,higherR)  - Filter(1,higherR) ;

if helpplot == 1
    figure(1234);subplot(1,2,1); 
    loglog(Fr,Filter(1,:),'color',[0.8 0.8 0.8],'linewidth',6);hold on;grid on
    loglog(Fr(~RTT & ~RTT2),Filter(1,~RTT & ~RTT2),'linewidth',2);
    loglog(Fr(RTT),Filter(1,RTT),'linewidth',2);
    loglog(Fr(RTT2),Filter(1,RTT2),'linewidth',2);
    %loglog(Fr,Filter(1,:),'k--');
    legend('all','Cardiff data','interm Freq patch','high freq patch','Location','southwest');
    xlim([1e11 1.5e14]);xlabel('Frequency  (Hz)');ylabel('S21');ylim([1e-6 2]);
    subplot(1,2,2);
    loglog(Fr,Filter(1,:));hold on;ylim([0 1])
    loglog(Fr,Filter(2,:));
    loglog(Fr,Filter(3,:));
    loglog(Fr,Filter(1,:)+Filter(3,:)+Filter(2,:),'k');
    legend('S21','S11','Absorption','Sum','Location','southwest');grid on;axis tight;ylim([1e-4 2])
    MakeGoodFigure(16,5,12,'filtercheck/FilterExample')
end

end

function [F,S21] = ReadaFilter(FN,wn)
%FN= asci tesxt full filename, 2 cols F(Hz) or wavenumber cm-1, S21 tab delimited
%wn = 1 then fisrt col is converted to F(Hz) from wavenumber
h = 6.6262e-34;		% J.s
c = 2.9998e8;		% m/s
k = 1.3806e-23;		% J/K

bba=dlmread(FN);%BPF Sample Holder
F=bba(:,1)';  %
S21=bba(:,2)';%patcghing 0 transmission at lowest F as this is a full BPF stack
if wn == 1
    F = F*100*c;% from cm-1 m-1 (*100) to Hz (*c)
end
end
