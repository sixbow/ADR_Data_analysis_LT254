function NoiseFilePatch
%patching files
%Reading FFT file
%Reading S21 file (GHz,Re,Im)
%Reading S21 file (GHz,dB,rad)
%Reading time-domain file
close all;     %close all open plots to remove clutter.
clc
addpath([pwd,filesep,'..',filesep,'subroutines']);

% specify where
MainDir =  '\\MARS\kid\KIDonSun\experiments\Entropy ADR\LT165_W2_Chip10\Noise_vs_T';
DirID = 'FFT';
SubDir = 'Power';
NewDir = [MainDir filesep 'Noise' filesep '2D'];%should not conatin 'FFT'
if ~isfolder(NewDir)
    mkdir(NewDir);
end

dirs=dir([MainDir filesep '*' DirID '*']);

nn=1;
dirnum=zeros(1,length(dirs));
for n=1:length(dirs)
    if dirs(n).isdir==1
        nn=nn+1;
    end
end

noDirs=nn-1;
dirnum(noDirs+1:end)=[];
dirnum=sort(dirnum);
if isempty(dirnum)
    error('no data in specified directory or directore not found');
end
clc

%now we can read in all files if need be, for all folders
%% FFT file
ncol = 3;
Filen_ID = 'dBm__FFT.dat'; %exact end of the filename for name construiction
dopatch(Filen_ID, MainDir, SubDir, NewDir, noDirs, dirs, ncol)
%S21 file
ncol = 3;
Filen_ID = 'dBm__S21.dat'; %exact end of the filename for name construiction
dopatch(Filen_ID, MainDir, SubDir, NewDir, noDirs, dirs, ncol)
% S21dB file
ncol = 3;
Filen_ID = 'dBm__S21dB.dat'; %exact end of the filename for name construiction
dopatch(Filen_ID, MainDir, SubDir, NewDir, noDirs, dirs, ncol)
%td file
ncol = 2;
Filen_ID = 'dBm__td.dat'; %exact end of the filename for name construiction
dopatch(Filen_ID, MainDir, SubDir, NewDir, noDirs, dirs, ncol)
end

function dopatch(Filen_ID, MainDir, SubDir, NewDir, noDirs, dirs, ncol)
% get to the first dir, and get all files and info

dirtolook = [MainDir filesep dirs(1).name filesep SubDir filesep];
disp(['looking in ' dirtolook]);
RawFFTfiles = dir([dirtolook '*' Filen_ID]);
for nn=1:length(RawFFTfiles)
    if isempty(RawFFTfiles(nn).name)==1
        error('file empty');
    end
    temp=cell2mat(textscan(RawFFTfiles(nn).name,'%*3s %f %*c %f %*s'));
    Fns(nn,1)=temp(1);
    Fns(nn,2)=temp(2);
end


for nn = 1 : length(Fns) %filenames
    for mm = 1 : noDirs %all dirs
        if mm == 1 %start dir
            %FFT file
            dirtolook = [MainDir filesep dirs(mm).name filesep SubDir];
            dirs(mm).name
            FFTfn = ['KID' num2str(Fns(nn,1)) '_' num2str(Fns(nn,2)) Filen_ID ];
            [Data,Temperature,Power,header] = import_data([dirtolook filesep FFTfn]);
            %write new file
            fid = fopen([NewDir filesep FFTfn],'wt');
            for m = 1:length(header)
                fprintf(fid, '%s\n',header{m});
            end
            fprintf(fid, 'Temperature in K:%.7f\n',Temperature);
            fprintf(fid, 'I	Q	dt= 0.000020\n');%this MUST be here, only td file relevant
            if ncol == 3
                fprintf(fid, '%.9f %.9f %.9f\n',Data{1}');
            elseif ncol == 2
                fprintf(fid, '%.9f %.9f\n',Data{1}');
            else
                error('ncol wrong')
            end
            clear Data
        else
            %Append only
            dirtolook = [MainDir filesep dirs(mm).name filesep SubDir ];
            dirs(mm).name
            if isfile([dirtolook filesep FFTfn])
                [Data,Temperature,Power] = import_data([dirtolook filesep FFTfn]);
                fprintf(fid,'\nTemperature in K:%.7f\n',Temperature);
                fprintf(fid,'Hz	Phase noise	Amp noise\n');
                if ncol == 3
                    fprintf(fid, '%.9f %.9f %.9f\n',Data{1}');
                elseif ncol == 2
                    fprintf(fid, '%.9f %.9f\n',Data{1}');
                else
                    error('ncol wrong')
                end
                clear Data
            end
        end
        
    end %end filename loop
    fclose(fid);
end %end dir loop
end
