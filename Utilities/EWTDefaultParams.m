function params = EWTDefaultParams()

% Wavelet type (littlewood-paley,shannon,gabor1 (rays=plain Gaussians),
% gabor2 (rays=half Gaussian/half cst)) for 1D transform.
params.wavname = 'littlewood-paley';


% Perform the detection on the log spectrum instead the spectrum
params.log=0;

params.SamplingRate = -1; %put -1 if you don't know the sampling rate

% Choose the wanted preprocessing for the scales (none,plaw,poly,morpho,tophat)
params.globtrend = 'none';
params.degree=5; % degree for the polynomial interpolation

% Choose the wanted regularization (none,gaussian,average,closing)
params.reg = 'gaussian';
params.lengthFilter = 10;
params.sigmaFilter = 1.5;

% Choose the wanted detection method for the scales (locmax,locmaxmin,ftc,scalespace)
params.detect = 'scalespace';
params.typeDetect='otsu'; %for scalespace:otsu,halfnormal,empiricallaw,mean,kmeans

params.N = 4; % maximum number of band for the locmaxmin method
params.completion=0; % complete (1) or not (0) the number of boundaries to reach the wanted number.


%% Used 1D Wavelet
%choose wavelet (littlewood-paley,shannon,gabor1 (rays=plain Gaussians),
% gabor2 (rays=half Gaussian/half cst))
params.wavname = 'littlewood-paley';


%% CURVELET PARAMETERS
% Type of curvelet transform (1=scale and radius independent, 2=angles per
% scales, 3=scales per angles)
params.option=3;

% Choose the wanted regularization for the angles
% (none,gaussian,average,closing)
params.curvreg = 'none';
params.curvlengthFilter = 10;
params.curvsigmaFilter = 1.5;

% Choose the wanted global trend removing for the angles (none,plaw,poly,morpho,tophat)
params.curvpreproc='none';
params.curvdegree=4;

% Choose the wanted detection method for the angles (locmax,locmaxmin)
params.curvmethod='scalespace';
params.curvN=2;

%% WATERSHED EWT PARAMETERS
%scale-space parameters
params.t = 0.8;
params.n = 6;
params.niter = 4;
params.nlevel = 100;    %for nlevel option
    
%transition width parameter
params.tau = .1;
    
%maxima post-processing parameters
params.includeCenter = 0;
params.edges = 0;

%transform type (0=real/1=complex)
params.complex = 0;
