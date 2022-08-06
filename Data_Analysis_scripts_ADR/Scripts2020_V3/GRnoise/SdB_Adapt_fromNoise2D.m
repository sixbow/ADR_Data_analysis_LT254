% Data_analysis LT254
%
%DATE: October 24, 2012 - Aug 20 2019
%Adapted by: Sietse de Boer
%AUTHOR of source NOISEanalysis_2DV1: Reinier Janssen, Jochem Baselmans
%
%==========================================================================
close all
clear all
global FITF0
%==========================================================================
% Inputs
%==========================================================================
ChipInfo.path = ['../..' filesep]; %root path where data is, one higher than the scripts
%FFTsubdir = ['Noise 120mK' filesep 'FFT' filesep 'Power'];  120mK    %
FFTsubdir = ['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D'];     %

%\\FFT\2D
ChipInfo.IndexPref = 1;     %Index of the power that is used as reference for Popt finding. (should not matter here)
numlowPremoved = 0;         %for noise plots(P) - removes lowest power for less clutter
SaveStuff = 1;              %0 to not save files (faster) 1 to save
Tsystem = 3;                %estimated system noise, 8K oid
BWrange = 1;                %Width aroung Fres of the reference power that is considered to find Popt. Default 1
FITF0 = 0;                    %Switch used to determine which method is used for determination of F0.
MaxnT = 18;     %max. #temperatures to be read.

%Please see the FitS21main5 routine for details. FITF0 = 0 is recommended due to the presence of overdriven
%resonators in Power sweeps (it will only update Q, not F0. F0 mot very important here).

%==========================================================================
% fixed settings
ReadPoptfile = 1;             % Popt must be obtained by doing NOISEanalysis_Pdep first!

%==========================================================================
% Start by setting some bla
%==========================================================================
format('long','g');     %Set display format of numbers to 7 digits
addpath([pwd,filesep,'..',filesep,'subroutines']);
close all;              %close all open plots to remove clutter.

%==========================================================================
% read Poptfile.
%==========================================================================
Noisepath = [ChipInfo.path,FFTsubdir,filesep]; %Path containing the raw noise data
disp(['Searchpath is: ' Noisepath])
if ReadPoptfile==1
    PoptFile = [Noisepath,'Popt.csv'];
    fid = fopen(PoptFile);
    if fid == -1 %no Popt file
        error([Noisepath,'Popt.csv not found, please do P analysis first'])
    else %There is a Popt file. Read it and use it.
        fclose(fid);
        fprintf('Found Popt.csv, this will be used to indicate Popt.\n')
        [~,PoptData] = ReadSRONcsvV2(PoptFile,'',0);    %PoptData will be Nx3, the cols are [T,KIDid,|Popt [dBm]|]
    end
end
clear PoptFile FFTsubdir
%==========================================================================
%Search the noisepath for all FFT files and filter out all files that are
%significantly longer or shorter than the mean.
%==========================================================================

RawFFTfiles = dir([Noisepath,'KID*FFT*.dat']);
if isempty(RawFFTfiles)
    error(['No data available in' Noisepath]);
end

%Determine the KID numbers (IDs) from the names of the good files.
KIDnumbers = zeros(1,length(RawFFTfiles));
NOISE(length(RawFFTfiles)).KIDnumber = ''; %pre allocate
for p=1:length(RawFFTfiles) %Loop over all FFT files
    KIDnumbers(p) = cell2mat(textscan(RawFFTfiles(p).name,'%*3s %f %*s'));
    NOISE(p).KIDnumber = KIDnumbers(p);
end
KIDnumbers = unique(KIDnumbers); %Determine all unique KIDs.

%Print to screen some of the reading information.
fprintf('Search for FFT data performed in:\n')
disp(Noisepath)
fprintf(['Inside the path the following KIDs are available: ',num2str(KIDnumbers),'\n'])
fprintf(['Inside this path a total number of ',num2str(length(RawFFTfiles)),' files were found.\n'])

%==========================================================================
% Read in all the data files
%==========================================================================
nTemp = zeros(length(RawFFTfiles),1);
for p=1:length(RawFFTfiles) %LOOP OVER ALL FILES (aka KID-P-combinations)
    LocEndRoot = strfind(RawFFTfiles(p).name,'FFT.dat');
    if isempty(LocEndRoot)
        fprintf('ERROR Noise Analysis: Cannot find FFT.dat in the name of the noise file.\n')
        fprintf('Most likely text has been placed between FFT and .dat extension.\n')
        error('Cannot create RootName for file detection.\n')
    end
    
    RootName = [Noisepath,RawFFTfiles(p).name(1:LocEndRoot-1)];
    NOISE(p).filename = RootName;
    clear LocEndRoot RootName
    %======================================================================
    %One by one read in the data files
    %======================================================================
    
    %======================================================================
    %Reading FFT file
    FFTfile = [NOISE(p).filename,'FFT.dat'];
    [Data,Temperature,Power,FFTheader] = import_data(FFTfile);%Data = cell array
    % limit T range
    if p == 1
        disp('Temperatures:')
        disp(Temperature);
    end
    % limit T range to max
    if length(Temperature) > MaxnT
        Temperature(MaxnT+1:end)=[];
    end
    MaxnT = length(Temperature);
    NOISE(p).Temperature = Temperature;%
    Data(MaxnT+1:end)=[];
    NOISE(p).FFTnoise = Data;
    clear Data;
    NOISE(p).ReadPower = -1*Power;
    nTemp(p) = MaxnT;
    
    %Find in FFTheader the frequency of the tone used to read the noise.
    for hl = 1:length(FFTheader)
        Freadindex = cell2mat(strfind(FFTheader(hl,1),'F used')); %Try to find F used in this line
        ColonIndex = cell2mat(strfind(FFTheader(hl,1),':')); %Check for : in this line
        if isempty(Freadindex)
            %Not in this line
        else
            %If it is in the line, find the first : after F used.
            Freadstart = find(ColonIndex>Freadindex,1);
            Freadstart = ColonIndex(Freadstart);
            Fread = cell2mat(textscan(FFTheader{hl,1}(Freadstart+1:end),'%f'));
        end
    end
    NOISE(p).Fread = Fread;
    clear ColonIndex Fread Freadindex Freadstart FFTfile hl FFTheader
    
    %======================================================================
    %Reading S21 file (GHz,Re,Im)
    S21file = [NOISE(p).filename,'S21.dat'];
    [Data,~,Power,~] = import_data(S21file);
    %Check for equality of temperature and power
    if abs(-1*Power-NOISE(p).ReadPower) > 0.5
        disp('Warning: Read Power does not match between FFT and S21 file.')
    end
    %Copy data for struct
    Data(MaxnT+1:end)=[];
    NOISE(p).S21_IQplane = Data;
    clear S21file Data
    %======================================================================
    %Reading S21 file (GHz,dB,rad)
    S21DBfile = [NOISE(p).filename,'S21dB.dat'];
    [Data,~,Power] = import_data(S21DBfile);
    %Check for equality of temperature and power
    if abs(-1*Power-NOISE(p).ReadPower(1)) > 0.5
        disp('Warning: Read Power does not match between FFT and S21dB file.')
    end
    %Copy data for struct
    Data(MaxnT+1:end)=[];
    NOISE(p).S21_MPplane = Data;
    clear Data S21DBfile
    %======================================================================
    %Reading time-domain file
    TDfile = [NOISE(p).filename,'td.dat'];
    [Data,~,Power,TDHeader] = import_data(TDfile);
    %Check for equality of temperature and power
    if abs(-1*Power-NOISE(p).ReadPower(1)) > 0.5
        disp('Warning: Read Power does not match between FFT and TD file.')
    end
    %Extract the timestep from the first subheader.
    fid=fopen(TDfile);
    dt = cell2mat(textscan(fid,'%*s%*s%*s%f',1,'headerlines',length(TDHeader)+1));
    fclose(fid);
    %Copy data for struct
    for nT=1:MaxnT
        CompleteData = zeros(size(Data{1},1),size(Data{1},2)+1);
        CompleteData(:,1) = dt*(1:size(Data{nT},1))';
        CompleteData(:,2:3) = Data{nT}(:,:);
        NOISE(p).TDIQ{nT,1} = CompleteData;
    end
    
    disp(['Highest Read: ' num2str(NOISE(p).Temperature(end))]);
    clear Temperature Data dt TDfile CompleteData fid TDHeader
    
    %======================================================================
    %Analyse the S21 data measured in [mag,phase] space for each temperature
    %======================================================================
    for nT=1:length(NOISE(p).Temperature)
        %Normalize the S21 data measured in magnitude plane
        filtered=smooth(NOISE(p).S21_MPplane{nT}(:,2),3); %smoothing the |S21| data in dB space
        %normalise S21 in log space to the max(|S21|)
        NOISE(p).S21_MPplane{nT}(:,2)=NOISE(p).S21_MPplane{nT}(:,2)-max(filtered);
        %Convert dB to magnitude
        NOISE(p).S21_MPplane{nT}(:,2) = 10.^(NOISE(p).S21_MPplane{nT}(:,2)/20);
        %Perform fit to obtain resonator parameters. Note the FITF0 is
        %recommended for overdriven resonators.
        [fres,Q,S21min,FitResult] = FitS21main5(NOISE(p).S21_MPplane{nT}(:,1:3),FITF0);
        %Put the temporary storage variables into the NOISE struct
        NOISE(p).Ql(nT)=Q(1);
        NOISE(p).Qi(nT)=Q(2);
        NOISE(p).Qc(nT)=Q(3);
        NOISE(p).Fres(nT)=fres;
        NOISE(p).S21min(nT)=S21min;%in dB!
        NOISE(p).S21fit{nT,1} = FitResult;
        %Calculate the Internal Power
        NOISE(p).InternalPower(nT) = 10*log10((2/pi)*10.^(NOISE(p).ReadPower/10).*(NOISE(p).Ql(nT).^2./NOISE(p).Qc(nT)));
    end
    %Calculate the resonator ring time
    NOISE(p).TauRes = NOISE(p).Ql(nT)./(pi*NOISE(p).Fres);
    %Calculate the resonator bandwidth (used later in Popt determination)
    NOISE(p).Bandwidth = NOISE(p).Fres./NOISE(p).Ql;
    %NOISE(p).MeanNoise = zeros(MaxnT,3);
    %NOISE(p).MeanFreqNoise = zeros(MaxnT,4);
    
    for nT=1:length(NOISE(p).Temperature)
        %Calculate frequency noise
        %Normalized Frequency Noise (Sf/F^2) [1/Hz] (defined by J. Gao) - setup
        %corrected
        tempnoise = (10.^(NOISE(p).FFTnoise{nT}(:,2)/10))*(1/(4*NOISE(p).Ql(nT))^2);
        Setup_level = mean(tempnoise(end-10:end-5));
        NOISE(p).FFTnoise{nT}(:,4) = tempnoise - Setup_level;
        %NOISE(p).FFTnoise{1}(NOISE(p).FFTnoise{1}(:,4)<0,4) = NaN;
        clear tempnoise Setup_level fres Q S21min FitResult
        
        %Frequency Noise (defined by B. Mazin) not used
        %NOISE(p).FFTnoise{nT}(:,5) = NOISE(p).FFTnoise{nT}(:,4)*(NOISE(p).Fres(nT)*1e9)^2;
        
        %Phase noise WRT complex plane
        S21min_a=10^(( NOISE(p).S21min(nT) )/20);%from dB to magnitude, not stored
        R=20*log10((1-S21min_a)/2);%correction due to circl radius vs complex plane radius in dB
        NOISE(p).FFTnoise{nT}(:,6)=NOISE(p).FFTnoise{nT}(:,2)+R;
        
        %Amplitude noise WRT complex plane
        NOISE(p).FFTnoise{nT}(:,7)=NOISE(p).FFTnoise{nT}(:,3)+R;
        
        %         %Calculate Mean amplitude noise between 20 and 500 Hz (Required for
        %         %Popt determination)
        %         Frange = 20 <= NOISE(p).FFTnoise{nT}(:,1) & NOISE(p).FFTnoise{nT}(:,1) <= 500;
        %         NOISE(p).MeanNoiseforPopt(nT) = mean(NOISE(p).FFTnoise{nT}(Frange,3));
        
        %Calculate the mean frequency noise (Gao definition), with setup
        %contribution subtracted (was done already in defining Gao's freq
        %noise)
        if min(NOISE(p).FFTnoise{(nT)}(:,1)) < 1000 %data exists (no cosmiuc ray proble,)
            NOISE(p).MeanFreqNoise1Hz(nT) = ... % 1 Hz
                mean(interp1(NOISE(p).FFTnoise{(nT)}(:,1),NOISE(p).FFTnoise{(nT)}(:,4),(0.8:0.1:1.2),'linear',1));
            NOISE(p).MeanFreqNoise10Hz(nT) = ... % 10 Hz
                mean(interp1(NOISE(p).FFTnoise{(nT)}(:,1),NOISE(p).FFTnoise{(nT)}(:,4),(9:0.1:11),'linear',1));
            NOISE(p).MeanFreqNoise100Hz(nT) = ... % 100 Hz
                mean(interp1(NOISE(p).FFTnoise{(nT)}(:,1),NOISE(p).FFTnoise{(nT)}(:,4),(90:1:110),'linear',1));
            NOISE(p).MeanFreqNoise1000Hz(nT) = ... % 1 kHz
                mean(interp1(NOISE(p).FFTnoise{(nT)}(:,1),NOISE(p).FFTnoise{(nT)}(:,4),(900:10:1100),'linear',1));
        else
            NOISE(p).MeanFreqNoise1Hz(nT) = NaN;
            NOISE(p).MeanFreqNoise10Hz(nT) = NaN;
            NOISE(p).MeanFreqNoise100Hz(nT) = NaN;
            NOISE(p).MeanFreqNoise1000Hz(nT) = NaN;
        end
    end
    
    
end %END OF LOOP OVER ALL FILES (aka KID-Pread-combinations)
clear Power R S21min_a Data Frange RawFFTfiles FITF0

% Clean tempertures. Sometimes not all files have the same number of
% temperatures
if ~all(nTemp == nTemp(1))
    MaxnT = min(nTemp);
    for p=1:length(NOISE)
        if length(NOISE(p).Temperature) > MaxnT
            NOISE(p).Temperature(MaxnT+1:end) = [];
            NOISE(p).MeanFreqNoise1Hz(MaxnT+1:end) = [];
            NOISE(p).MeanFreqNoise10Hz(MaxnT+1:end) = [];
            NOISE(p).MeanFreqNoise100Hz(MaxnT+1:end) = [];
            NOISE(p).MeanFreqNoise1000Hz(MaxnT+1:end) = [];
            NOISE(p).Temperature(MaxnT+1:end) = [];
            NOISE(p).Bandwidth(MaxnT+1:end) = [];
            NOISE(p).TauRes(MaxnT+1:end) = [];
            NOISE(p).InternalPower(MaxnT+1:end) = [];
            NOISE(p).S21fit(MaxnT+1:end) = [];
            NOISE(p).S21min(MaxnT+1:end) = [];
            NOISE(p).Fres(MaxnT+1:end) = [];
            NOISE(p).Qc(MaxnT+1:end) = [];
            NOISE(p).Qi(MaxnT+1:end) = [];
            NOISE(p).Ql(MaxnT+1:end) = [];
            NOISE(p).S21fit(MaxnT+1:end) = [];
            NOISE(p).TDIQ(MaxnT+1:end) = [];
            NOISE(p).S21_MPplane(MaxnT+1:end) = [];
            NOISE(p).S21_IQplane(MaxnT+1:end) = [];
            NOISE(p).FFTnoise(MaxnT+1:end) = [];
        end
    end
end

%==========================================================================
% Now that all files have been read we read Popt fro csv file, 1 value for
% all temperatures
%==========================================================================
IndexPsort = cell(length(KIDnumbers),1);
% initialize cell array to store for each KID an array that contains:
% contains vector of length = number of readout powers (NP). This vector holds the indexes in NOISE
% that contain information about this resonator. The indexes are given in
% order of increasing power.
% IndexPsort{kidn} = array all indices for kidn in NOISE

NOISEParams = cell(length(KIDnumbers),1);IndexP_sub_opt = NOISEParams;IndexP_super_opt = NOISEParams;
IndexPopt = zeros(length(KIDnumbers),1);
for kidn=1:length(KIDnumbers) % LOOP OVER ALL UNIQUE KIDS,
    %======================================================================
    % create IndexPsort{kidn}: indices with Pread sorted,, valid for all T
    %======================================================================
    KIDlocations = find([NOISE(:).KIDnumber] == KIDnumbers(kidn));  %gets the indices in NOISE that contains data on this KID (i.e which files this KID is descrived in)
    [~,blaindices]=sort([NOISE(KIDlocations).ReadPower]);           %Sort on Pread and save the indexes.
    IndexPsort{kidn} = KIDlocations(blaindices);                    %Recalibrate the indexes to the NOISE indices. This vector can be used to get the sorted P read for this KID, valid for all temperaturesd
    clear KIDlocations blaindices
    %Check if the desired reference index is within the number of powers measured for this KID%
    if ChipInfo.IndexPref > length(IndexPsort{kidn})
        disp('WARNING NoiseAnalysis: IndexPref exceeds number of Powers. Assuming maximum Power to be reference.')
        ChipInfo.IndexPref = length(IndexPsort{kidn});
    end
    
    
    
    %local Shorthand variable for the reference power index (valid locally for this kid as we loop over KIDs here)%
    IndexPref = IndexPsort{kidn}(ChipInfo.IndexPref);
    
    %==================================================================
    % DETERMINE Popt
    %==================================================================
    %Determine which powers have to big a resonance frequency shift
    if ReadPoptfile == 0
        error('Popt.csv not found; run P analysis on 2D data first')
    elseif ReadPoptfile==1
        %Put some user info to screen. ReadPower(1) refers to lowest T
        disp(['Reading optimum power for KID ',num2str(KIDnumbers(kidn))]);
        %
        Poptdata_Tindices= PoptData(:,1)>=0.9*NOISE(IndexPref).Temperature(1,1) & PoptData(:,1)<=1.1*NOISE(IndexPref).Temperature(1,1);%Note MaxnT = length(NOISE(IndexPref).Temperature);
        Poptdata_thisKID_thisT=PoptData(Poptdata_Tindices,2)==NOISE(IndexPref).KIDnumber;
        bla=find(Poptdata_Tindices);%this is an array of indiced in PoptData that are allowed wrt temperature
        Popt_fromfile=PoptData(bla(Poptdata_thisKID_thisT),3);
        tempPoptIndex=find([NOISE(IndexPsort{kidn}).ReadPower]==Popt_fromfile);
    else
        error('ReadPoptfile must be 0 or 1')
    end
    %==================================================================
    % New indices for NOISE struct: One for Popt, one for all powers below Popt. Same for all T%
    %==================================================================
    % ******************* IndexPopt *******************
    IndexPopt(kidn) = IndexPsort{kidn}(tempPoptIndex); %Index of Popt in full NOise struct (for each T)  IndexPopt(kidn)
    
    % find the indices above and below in IndexPsort{kidn}
    temp_sub_opt    = ([NOISE(IndexPsort{kidn}).ReadPower] <= NOISE(IndexPopt(kidn)).ReadPower); %index in Indexpsort
    temp_super_opt  = ([NOISE(IndexPsort{kidn}).ReadPower] > NOISE(IndexPopt(kidn)).ReadPower);
    
    % ******************* indices above and below popt *******************
    IndexP_sub_opt{kidn}    = IndexPsort{kidn}(temp_sub_opt);%Index of all powers below Popt in full NOise struct (for each T)  IndexPopt(kidn)
    IndexP_super_opt{kidn}  = IndexPsort{kidn}(temp_super_opt);
    
    % create Noise paraemeters struct to allow eacy KID to KID plots
    NOISEParams{kidn}.Popt = NOISE(IndexPopt(kidn)).ReadPower;
    NOISEParams{kidn}.Popt_int = NOISE(IndexPopt(kidn)).InternalPower;
    clear FresPref AllowedPindices PoptIndex tempPoptIndex Poptdata_Tindices Poptdata_thisKID_thisT bla Popt_fromfile temp_sub_opt temp_super_opt
    
    
    %==================================================================
    % Figure: T dependencies for every power. Get KID number and Power
    % right sing Prr
    %==================================================================
    Pcolors = colormapJetJB(length(IndexPsort{kidn}));
    Tcolors = colormapJetJB(MaxnT);
    
    for Prr = 1:length(IndexP_sub_opt{kidn})
        
        Pr = IndexP_sub_opt{kidn}(Prr);
        figure(1000*NOISE(Pr).KIDnumber+Prr)
        clf
        
        clear TempLegend
        TempLegend = cell(length(NOISE(Pr).Temperature)+1,1);
        TempLegend{1} = 'P_{opt}';
        
        %Resonance circle as a function of T (incl noise blobs)
        subplot(2,3,1)
        hold on
        for nT=1:length(NOISE(Pr).Temperature)
            plot(NOISE(Pr).S21_IQplane{nT}(:,2) , NOISE(Pr).S21_IQplane{nT}(:,3),'-','color',Tcolors(nT,:),'LineWidth',1)
            TempLegend{nT+1} = num2str(NOISE(Pr).Temperature(nT),'%.3f');
        end
        for nT=1:length(NOISE(Pr).Temperature) %Second Loop to get the legend correct.
            plot(NOISE(Pr).TDIQ{nT}(:,2),NOISE(Pr).TDIQ{nT}(:,3),'.','color',Tcolors(nT,:),'MarkerSize',6)
        end
        legend(TempLegend(2:end))
        xlabel('Re');ylabel('Im')
        title(['KID ',num2str(NOISE(Pr).KIDnumber(1),'%.0f'),' @P_{read}=',num2str(NOISE(Pr).ReadPower,'%.0f'),' dBm'])
        box on;grid on;
        hold off
        
        %Resonance Dip as a function of T, Incl reference lines around reference power%
        subplot(2,3,2)
        hold on
        for nT=1:length(NOISE(Pr).Temperature)
            plot(NOISE(Pr).S21_MPplane{nT}(:,1),20*log10(NOISE(Pr).S21_MPplane{nT}(:,2)),'-','color',Tcolors(nT,:),'LineWidth',1)
        end
        xlabel('F [GHz]');ylabel('|S21| [dB]')
        title(['KID ',num2str(NOISE(Pr).KIDnumber(1),'%.0f'),' P_{opt}=',num2str(NOISE(Pr).ReadPower),' dBm'])
        axis tight;
        box on;grid on;
        hold off
        
        %Time Domain Trace as a function of time
        subplot(2,3,3)
        semilogy(NOISE(Pr).Temperature,NOISE(Pr).Qi,'ro','MarkerSize',6,'MarkerFaceColor','r');%At Popt
        xlabel('T  (K)');ylabel('Q_{i} (dBm)')
        box on;grid on;
        hold off
        
        %Frequency Noise
        subplot(2,3,4)
        warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart');
        toplot = NOISE(Pr).FFTnoise{1}(:,4) > 0;
        semilogx(NOISE(Pr).FFTnoise{1}(toplot,1),10*log10(NOISE(Pr).FFTnoise{1}(toplot,4)),...
            '-','color','k','LineWidth',2)
        hold on
        for nT=1:length(NOISE(Pr).Temperature)
            toplot = NOISE(Pr).FFTnoise{nT}(:,4) > 0;
            semilogx(NOISE(Pr).FFTnoise{nT}(toplot,1),10*log10(NOISE(Pr).FFTnoise{nT}(toplot,4)),...
                '-','color',Tcolors(nT,:),'LineWidth',1)
        end
        xlabel('F [Hz]');ylabel('S_F/F^2 [dBc/Hz]')
        xlim([0.5,1e5]);grid on;ylim([-220,-140])
        hold off
        
        %Phase Noise (and Amp Noise)
        subplot(2,3,[5,6])
        semilogx(NOISE(Pr).FFTnoise{1}(:,1),NOISE(Pr).FFTnoise{1}(:,3),...
            '-','color','k','LineWidth',3)
        hold on
        semilogx(NOISE(Pr).FFTnoise{1}(:,1),NOISE(Pr).FFTnoise{1}(:,2),...
            '-','color','k','LineWidth',4)
        minnp = zeros(length(NOISE(Pr).Temperature),1);maxnp = minnp;
        for nT=1:length(NOISE(Pr).Temperature)
            % finding range
            maxnp(nT) = max(NOISE(Pr).FFTnoise{nT}(NOISE(Pr).FFTnoise{nT}(:,1) > 10,2));
            minnp(nT) = mean(NOISE(Pr).FFTnoise{nT}(20:end-5,3));
            semilogx(NOISE(Pr).FFTnoise{nT}(:,1),NOISE(Pr).FFTnoise{nT}(:,3),...
                '-','color',Tcolors(nT,:),'LineWidth',1)
            semilogx(NOISE(Pr).FFTnoise{nT}(:,1),NOISE(Pr).FFTnoise{nT}(:,2),...
                '-','color',Tcolors(nT,:),'LineWidth',2)
        end
        grid on;
        xlabel('F [Hz]');ylabel('S_x [dBc/Hz]')
        legend('S_A','S_{\theta}')
        xlim([3,0.3e6]);ylim([10*floor(min(0.1*minnp))-5,10*ceil(max(0.1*maxnp))+5]);
        hold off
        
        %SAVE the figure
        Figfile=[Noisepath,'KID',num2str(NOISE(Pr).KIDnumber,'%.0f'),'_Pr_',...
            num2str(-1*NOISE(Pr).ReadPower,'%.3g'),'dBm_NOISE_vsT'];
        if SaveStuff == 1
            MakeGoodFigure(15,15,14,Figfile)
        else
            MakeGoodFigure(15,15,14)
        end
    end %power loop
    close all;
end %END OF LOOP OVER ALL UNIQUE KIDS


rmpath([pwd,filesep,'..',filesep,'subroutines']);
if SaveStuff == 1
    save([Noisepath,'Noise_2D.mat'], '-v7.3')
end
