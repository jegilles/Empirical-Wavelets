%% Script to test the Empirical Wavelet Transform
% Based on the papers:
% J. Gilles, "Empirical Wavelet Transform", IEEE
% Trans. on Signal Processing, 2013
% J. Gilles, G. Tran, S. Osher, "2D Empirical transforms. Wavelets, 
% Ridgelets and Curvelets Revisited", SIAM Journal on Imaging Sciences,
% 2014
% J. Gilles, K. Heal, "A parameterless scale-space approach to find 
% meaningful modes in histograms - Application to image and spectrum 
% segmentation", submitted 2014.
%
% Don't hesitate to modify the parameters and try
% your own signals!
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2014
% Version: 2.0
clear all
close all

%% User setup

% Choose the signal you want to analyze
% (sig1,sig2,sig3,sig4=ECG,sig5=seismic,sig6=EEG,sig7=compplex signal1)
signal = 'sig3';
params.SamplingRate = -1; %put -1 if you don't know the sampling rate
%params.SamplingRate = 4000; %put -1 if you don't know the sampling rate
%channel = 50; %for EEG only

%choose wavelet (littlewood-paley,shannon,gabor1 (rays=plain Gaussians),
% gabor2 (rays=half Gaussian/half cst))
params.wavname = 'littlewood-paley';

% Choose the wanted global trend removal (none,plaw,poly,morpho,tophat,opening)
params.globtrend = 'none';
params.degree=6; % degree for the polynomial interpolation

% Choose the wanted regularization (none,gaussian,average,closing)
params.reg = 'none';
params.lengthFilter = 10;
params.sigmaFilter = 1.5;

% Choose the wanted detection method (locmax,locmaxmin,
% adaptive,adaptivereg,scalespace)
params.detect = 'scalespace';
params.typeDetect='otsu'; %for scalespace:otsu,halfnormal,empiricallaw,mean,kmeans
params.N = 4; % maximum number of bands
params.completion = 0; % choose if you want to force to have params.N modes
                       % in case the algorithm found less ones (0 or 1)
%params.InitBounds = [4 8 13 30];
params.InitBounds = [2 25];

% Perform the detection on the log spectrum instead the spectrum
params.log=0;

% Choose the results you want to display (Show=1, Not Show=0)
Bound=0;   % Display the detected boundaries on the spectrum
Comp=1;    % Display the EWT components
Rec=0;     % Display the reconstructed signal
TFplane=0; % Display the time-frequency plane (by using the Hilbert 
           % transform). You can decrease the frequency resolution by
           % changing the subresf variable below. (WORKS ONLY FOR REAL
           % SIGNALS
Demd=0;    % Display the Hilbert-Huang transform (WORKS ONLY FOR REAL 
           % SIGNALS AND YOU NEED TO HAVE FLANDRIN'S EMD TOOLBOX)
           
subresf=1;

InitBounds = params.InitBounds;
           
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
        load('eeg.mat')
        t=0:length(f)-1;
    case 'sig7'
        load('csig1.mat')
        t=0:length(f)-1;
end


%% We perform the empirical transform and its inverse
% compute the EWT (and get the corresponding filter bank and list of 
% boundaries)
[ewt,mfb,boundaries]=EWT1D(f,params);


%% Show the results

if Bound==1 %Show the boundaries on the spectrum
    div=1;
    if (strcmp(params.detect,'adaptive')||strcmp(params.detect,'adaptivereg'))
        Show_EWT_Boundaries(f,boundaries,div,params.SamplingRate,InitBounds);
    else
        Show_EWT_Boundaries(f,boundaries,div,params.SamplingRate);
    end
end

if Comp==1 %Show the EWT components and the reconstructed signal
        Show_EWT(ewt);
end

if Rec==1
        %compute the reconstruction
        rec=iEWT1D(ewt,mfb,~isreal(f));
        if isreal(f)
            figure;
            subplot(2,1,1);plot(f);title('Original signal');
            subplot(2,1,2);plot(rec);title('Reconstructed signal');
            disp('Reconstruction error:');
            norm(f-rec,Inf)
        else
            figure;
            subplot(2,2,1);plot(real(f));title('Original signal - real part');
            subplot(2,2,2);plot(imag(f));title('Original signal - imaginary part');
            subplot(2,2,3);plot(real(rec));title('Reconstructed signal - real part');
            subplot(2,2,4);plot(imag(rec));title('Reconstructed signal - imaginary part');
            disp('Reconstruction error:');
            norm(f-rec,Inf)            
        end
end

%Plot the time-frequency plane by using the Hilbert transform for the real
%case
if TFplane==1
    if isreal(f)
        EWT_TF_Plan(ewt,boundaries,params.SamplingRate,f,[],[],subresf,[]);
    else
        disp('TF plane plot not available for this type of signal.');
    end
end

%% EMD comparison: if you have Patrick Flandrin's EMD toolbox you can 
%% perform the EMD and display the corresponding Time-Frequency plane
if Demd==1
    imf=emd(f);
    Disp_HHT(imf,t,f,1,1);
end