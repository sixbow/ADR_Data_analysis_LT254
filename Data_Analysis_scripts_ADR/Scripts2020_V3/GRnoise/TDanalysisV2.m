% TDanalysisV2
% reads td binary data, calulates and plots the PSD's
% saves the data into a new crossPSD struct struct and appends to the exiting .mat file
% takes 11" on SRON1825 with local data per vfile
% you can save 5.5" by not plotting
clear all;close all;clc
tic
%NB: sometimes files are empty; delete those!
%NB: will NOT process binary files already analyzed.
%================================================================================
% Input
%================================================================================
ChipInfo_path = ['..' filesep '..' ]; %root path where data is, one higher than the scripts
PT2Ddep = 2;           %=0 for Tdep analysis, =1 for Pdep analysis, =2 for 2D analysis, error otheriwse
if PT2Ddep == 0
    %Tdep
    FFTsubsubdir=[ 'FFT' filesep 'Temp'];                   %FFTsubdir = [filesep 'Noise_Powers_165mK' filesep 'FFT' filesep 'Power'];     %
    TDsubsubdir=['FFT' filesep 'TD_Temp'];
elseif PT2Ddep == 1%Pdep
    FFTsubsubdir=['Noise 120mK' filesep 'FFT' filesep 'Power'];                   %FFTsubdir = [filesep 'Noise_Powers_165mK' filesep 'FFT' filesep 'Power'];     %
    TDsubsubdir=['Noise 120mK'  filesep 'TD_Power'];
elseif PT2Ddep == 2
    %FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D_Popt'];                   %FFTsubdir = [filesep 'Noise_Powers_165mK' filesep 'FFT' filesep 'Power'];     %
    FFTsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'FFT' filesep '2D'];                   %FFTsubdir = [filesep 'Noise_Powers_165mK' filesep 'FFT' filesep 'Power'];     %
    
    TDsubsubdir=['Data_LT254_Sietse' filesep 'LT254_Sietse_Chip11' filesep 'Noise_vs_T' filesep 'TD_2D'];
else
    error('PT2Ddep not defined')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WARNING: Set to 1 only if sure you have no data from previous analysis runs!!!%
killpreviousrunsforsure = 0;          % disables checking if files are already there (i.e. CrossPSDNOISE killed)
maakdeplots             = 1;            %1 makes plots, 0 disables plotting (2x faster almost)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timemed = 64;       % duration of td data @ 50 ksample/sec; is varied so check the files, 40 is old default = 31MB. 128"= 100 mB.
TDoption.nrsigma = 5;                 % 8 peaks>TDoption.nrsigma*sigma are rejected (based on filtered TD-stream) %6 seems to be really minimum, 8 good value also for smaller peaks
TDoption.smoothtime = 0.8e-3;           % 1e-3quasiparticle lifetime guess used for fit (note it will be updated if there is a resultfile present)
average_med = 128;                       % = TDoption.average for medium speed data. = 64,  (not more for 40'data)
% should be a power of 2 AND nrpoints/average_med/2 should be integer.
% Average is also the number of pieces over which peak rejection is done, ie if TDoption.average=32 and a peak is found,
% that 1/32th piece is thrown away.
% scipt uses TDp[tion.average, is set to 32  for the fast data and to average_med for med data%
TDoption.ppd = 30;                    % 30 points per decade in frequency for logsmooth after calculating PSD,
% 30 (20) will give spectrum of 158 (105) points. 0 will give you the full spectrum, usually a million points
% 30 is used standard in the SRON labview
% NOTE: if you want to use fitranges from earlier processing, make sure ppd is the same!!!
TDoption.savex = 0;                   % 0 save filtered time domain data in csv-file (note some 90 Mb per file); Saves in data direcotory
TDoption.corroff = 1;                 % 1 = correct each piece of length/average with linear offset (to correct for very slow drifts, only relevant for 50 kHz files
average_fast = 32; %do not touch

% set by HW VI's in labview
freqmed = 50e3;              %standard is 50e3;
freqfast = 1e6;              %standard is 1e6;
timefast = 0.2;              %standard is 0.2;

if rem( freqmed * timemed / average_med , 1)~=0
    error('TDoption.average = wrong for med data')
end
%================================================================================
if PT2Ddep == 0
    %Tdep
    matfile = 'Noise_T.mat';
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','KIDnumbers');
    matfile2 = 'CrossPSDNoise_T';
    if isfile([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile2])
        load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile2],'CrossPSDNOISE');
    end
    
elseif PT2Ddep == 1%Pdep
    matfile = 'Noise_P.mat';
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','IndexP_sub_opt','KIDnumbers','IndexPopt');
    matfile2 = 'CrossPSDNoise_P';
    if isfile([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile2])
        load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile2],'CrossPSDNOISE');
    end
    
elseif PT2Ddep == 2
    matfile = 'Noise_2D.mat';
    load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile],'NOISE','IndexP_sub_opt','KIDnumbers','IndexPopt');
    matfile2 = 'CrossPSDNoise_2D';
    if isfile([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile2])
        load([ChipInfo_path,filesep,FFTsubsubdir,filesep,matfile2],'CrossPSDNOISE');
    end
end

addpath([pwd,filesep,'..',filesep,'subroutines']);
addpath([pwd,filesep,'..',filesep,'Time Domain analysis']);

ChipInfo.path = ChipInfo_path;clear ChipInfo_path;

TDpath = [ChipInfo.path,filesep,TDsubsubdir,filesep]; %end with filesep
disp(['TDpath: ',TDpath]);

if ~exist('CrossPSDNOISE') % pre allocate
    CrossPSDNOISE(length(NOISE)).CrossPSD = '';
    disp('CrossPSDNOISE created')
end
if killpreviousrunsforsure == 1
    clear CrossPSDNOISE
    CrossPSDNOISE(length(NOISE)).CrossPSD = '';
end
toc
for kidn=1:length(KIDnumbers) % LOOP OVER ALL UNIQUE KIDS,
    %construct filename
    if PT2Ddep == 0
        IndexPopt(kidn) = kidn;
        IndexP_sub_opt{kidn} = kidn;
    end
    ID = num2str(NOISE(IndexPopt(kidn)).KIDnumber);
    for p=1:length(IndexP_sub_opt{kidn})% over Power
        for nT=1:length(NOISE(kidn).Temperature) % over T
            % construct filename
%             %Sietse Insert to find popt.
%             % KiD: [1 2 3 4 5 6]
%             P_Opt_order = [-94 -90 -94 -87 -92 -87];
%             p = 1;
%             while ~(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower == P_Opt_order(kidn))
%                p = p + 1; %What this does is it finds the right p for optimal power from the P_Opt_order list.
%             end
%             % Question how can i find index p such that
%             % NOISE(IndexP_sub_opt{kidn}(p)).ReadPower is at P_opt
%             % p should depend on kidn and nothing else
%             %/Sietse
            Pr = num2str(-1*NOISE(IndexP_sub_opt{kidn}(p)).ReadPower);
            Tk = num2str(round(1e3*NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT)),'%.0f');
            fnmed = ['KID' ID '_' Pr 'dBm__TDmed_TmK' Tk '.bin' ];  %bin filename - 50ksample/sec
            fnfast = ['KID' ID '_' Pr 'dBm__TDfast_TmK' Tk '.bin' ];%bin filename 2Msample/sec
            figname = [TDpath 'KID' ID '_' Pr 'dBm__TD_TmK' Tk];    %figure name for export
            fignamedisp = ['KID' ID '_' Pr 'dBm__TD_TmK' Tk];
            %check if files are there
            fid = fopen([TDpath fnmed ]);fid2 = fopen([TDpath fnfast ]);
            if fid == -1 || fid2 == -1
                disp(['WARNING import_data: Cannot find ' fnmed  ])
            else
                disp(['Reading: ' fnmed])
                fclose(fid);fclose(fid2); %close here; will open again in bimfilefunction4.m
            end
            
            %reading data, make figure; figure name defined as binfilefunction4. %
            % will only do so if figure file with the same name does NLT
            % exist
            if length(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD) >= nT % cell exists
                if ~isempty(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT})
                    disp([fnmed 'already analyzed, not done again'])
                end
            elseif fid == -1 || fid2 == -1
                disp('File skipped, NaN inserted in CrossPSD')
                CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT} = [];
                CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,4) = NaN;
                CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,3) = NaN;
                CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2) = NaN;
                CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1) = NaN;
            else
                if maakdeplots == 1
                    figure(1e7*NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber + 1e4*abs(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower) + round(1000*NOISE(kidn).Temperature(nT)));
                end
                [fmfast,SRRrawfast,SPPrawfast,SPRrawfast]   = binfilefunction4([TDpath fnfast ],timefast,freqfast,average_fast,TDoption, maakdeplots);
                [fmmed,SRRrawmed,SPPrawmed,SPRrawmed,~,~, CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).nrsigmaused(nT)] = ...
                    binfilefunction4([TDpath fnmed ],timemed,freqmed,average_med,TDoption, maakdeplots);
                %correcting and patching spectral data med to fast
                fmmed(end-5:end) = [];SRRrawmed(end-5:end) = [];SPPrawmed(end-5:end) = [];SPRrawmed(end-5:end) = [];
                fmfast(1:39) = [];SRRrawfast(1:39) = [];SPPrawfast(1:39) = [];SPRrawfast(1:39) = [];
                fm = [fmmed fmfast];SRRraw = [SRRrawmed SRRrawfast];SPPraw = [SPPrawmed SPPrawfast];SPRraw = [SPRrawmed SPRrawfast];
                %Put data in CrossPSDNOISE struct=
                CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT} = [];
                CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,4) = SRRraw;
                CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,3) = SPPraw;
                CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2) = SPRraw;
                CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1) = fm;
                
                %===================================================================
                %Plotting (cntnd from binfilefunction4)
                %===================================================================
                if maakdeplots == 1
                    subplot(2,2,3)%plot R and theta labview
                    semilogx(NOISE(IndexP_sub_opt{kidn}(p)).FFTnoise{nT}(:,1),NOISE(IndexP_sub_opt{kidn}(p)).FFTnoise{nT}(:,2),...
                        'r-','LineWidth',2);hold on;
                    semilogx(fm,SPPraw,'b-','LineWidth',2);
                    %add R
                    semilogx(NOISE(IndexP_sub_opt{kidn}(p)).FFTnoise{nT}(:,1),NOISE(IndexP_sub_opt{kidn}(p)).FFTnoise{nT}(:,3),'r-');
                    semilogx(fm,SRRraw,'b-');
                    grid on;axis tight;xlim([10,0.3e6]);
                    xlabel('F [Hz]');ylabel('S_x [dBc/Hz]');
                    legend('\theta Labview','\theta Matlab','R Labview','R Matlab')
                    
                    subplot(2,2,4)
                    semilogx(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,1),-1*real(CrossPSDNOISE(IndexP_sub_opt{kidn}(p)).CrossPSD{nT}(:,2)),'b-');grid on;
                    xlabel('F [Hz]');ylabel('S_{cross} [dBc/Hz]');
                    grid on;axis tight;xlim([10,0.1e6]);
                    title(['KID ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).KIDnumber) ' @Pread = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).ReadPower)...
                        ' dBm, T = ' num2str(NOISE(IndexP_sub_opt{kidn}(p)).Temperature(nT))]);
                    MakeGoodFigure(8,6,8,figname,1); %save png
                    close(gcf)
                    % save in the loop!
                end
                if isfile([ChipInfo.path,filesep,FFTsubsubdir,filesep,matfile2])
                    save([ChipInfo.path,filesep,FFTsubsubdir,filesep,matfile2],'CrossPSDNOISE','-append');
                else
                    save([ChipInfo.path,filesep,FFTsubsubdir,filesep,matfile2],'CrossPSDNOISE');
                end
            end
        end
    end
    toc
end
