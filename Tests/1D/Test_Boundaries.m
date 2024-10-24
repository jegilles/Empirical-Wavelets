%% Test of different Fourier boundaries detection strategies
clear all
close all

%% User setup
% Choose the signal you want to analyze
% (sig1,sig2,sig3,sig4=ECG,sig5=seismic,sig6=EEG,lena,textures)
signal = 'sig1';
params.SamplingRate = -1; %put -1 if you don't know the sampling rate

% Set to 1 to perform the detection on the log spectrum instead the spectrum
params.log=0;

% Choose the wanted global trend removal (none,plaw,poly,morpho,tophat)
params.globtrend = 'none';
params.degree=6; % degree for the polynomial interpolation
params.sizeel=4;

% Choose the wanted regularization (none,gaussian,average,closing)
params.reg = 'none';
params.lengthFilter = 5;
params.sigmaFilter = 1;

% Choose the wanted detection method (locmax,locmaxmin,localmaxminf,
% adaptive,adaptivereg,scalespace)
params.detect = 'scalespace';
params.typeDetect='otsu'; %for scalespace:otsu,halfnormal,empiricallaw,mean,kmeans
params.N = -1; % maximum number of band for the locmaxmin method
%params.alpha=0.25; %0.2 for eeg
params.completion = 0;

params.InitBounds = [10 25 75 100];
InitBounds = params.InitBounds; % keep a backup for visualization


%% Load signals
switch lower(signal)
    case 'sig1'
        load('sig1.mat');
        t=0:1/length(f):1-1/length(f);
    case 'sig2'
        load('sig2.mat');
        t=0:1/length(f):1-1/length(f);
    case 'sig3'
        load('sig3.mat');
        t=0:1/length(f):1-1/length(f);
    case 'sig4'
        load('sig4.mat');
        t=0:length(f)-1;
    case 'sig5'
        load('seismic.mat');
        f=f(10000:20000); %sub portion of the signal used in the paper
        t=0:length(f)-1;
    case 'sig6'
        load('eeg.mat') %directly loads f
        t=0:length(f)-1;
    case 'lena'
        load lena
        fftim=fft(f');
        ff=abs(sum(abs(fftim),2)/size(fftim,2));
        f=ifft(ff);
    case 'textures'
        load('texture.mat');
        fftim=fft(f');
        ff=abs(sum(abs(fftim),2)/size(fftim,2)); 
        f=ifft(ff);
end

if (~strcmp(signal,'lena')) && (~strcmp(signal,'textures')) && (~strcmp(signal,'eeg'))
    % We compute the Fourier transform of f
    ff=abs(fft(f));
end

%% Perform the detection and plot the detected boundaries
[boundaries,u] = EWT_Boundaries_Detect(ff(1:round(length(ff)/2)),params);
boundaries = boundaries*2*pi/length(ff);

div=1;
if (strcmp(params.detect,'adaptive')||strcmp(params.detect,'adaptivereg'))
    Show_EWT_Boundaries(f,boundaries,div,params.SamplingRate,params.log,InitBounds,u);
else
    Show_EWT_Boundaries(f,boundaries,div,params.SamplingRate,params.log);
end
