%% This script permits to test the 2D Curvelet EWT
% It generates all the results given in the paper
% J. Gilles, G. Tran, S. Osher, "2D Empirical tranforms. Wavelets, 
% Ridgelets and Curvelets Revisited" submitted at SIAM Journal on
% Imaging Sciences. 2013
%
% Don't hesitate to modify the parameters and try
% your own images!
% Warning: the image must be squared
%
% Author: Jerome Gilles - Giang Tran
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
clear

%% User setup
% Choose the image you want to analyze (texture,lena,barb)
signal = 'texture';

% Type of curvelet transform (1=scale and radius independent, 2=angles per
% scales, 3=scales per angles)
params.option=1;

% Choose the wanted preprocessing for the scales (none,plaw,poly,morpho,tophat)
params.globtrend = 'none';
params.degree=5; % degree for the polynomial interpolation

% Choose the wanted regularization (none,gaussian,average,closing)
params.reg = 'none';
params.lengthFilter = 10;
params.sigmaFilter = 1.5;

% Choose the wanted preprocessing for the angles (none,plaw,poly,morpho,tophat)
params.curvpreproc='none';
params.curvdegree=4;

% Choose the wanted regularization for the angles
% (none,gaussian,average,closing)
params.curvreg = 'none';
params.curvlengthFilter = 10;
params.curvsigmaFilter = 1.5;

% Choose the wanted detection method for the scales (locmax,locmaxmin,ftc,scalespace)
params.detect = 'scalespace';
params.typeDetect='otsu'; %for scalespace:otsu,halfnormal,empiricallaw,mean,kmeans

params.N = 4; % maximum number of band for the locmaxmin method
params.completion=0;

% Choose the wanted detection method for the angles (locmax,locmaxmin)
params.curvmethod='scalespace';
params.curvN=2;

% Perform the detection on the log spectrum instead the spectrum
params.log=0;

% Choose the results you want to display (Show=1, Not Show=0)
Bound=0;   % Display the detected boundaries on the spectrum
Comp=1;    % Display the EWT components
Rec=1;     % Display the reconstructed signal
           
switch lower(signal)
    case 'texture'
        load('texture.mat');
        %f=noisy;
    case 'lena'
        load('lena.mat');
    case 'barb'
        load('barb.mat');
end


%% We perform the 2D Littlewood-Paley EWT
[ewtc,mfb,Bw,Bt]=EWT2D_Curvelet(f,params);

%% Show results
if Comp==1 %Show curvelet components
   Show_EWT2D_Curvelet(ewtc,params.option);
end

if Rec==1 %Show the reconstructed image
   figure;subplot(1,2,1);imshow(f,[]);title('Original image');
   rec=iEWT2D_Curvelet(ewtc,mfb);
   disp('Reconstruction error:');
   norm(f(:)-rec(:),Inf)
   subplot(1,2,2);imshow(rec,[]);title('Reconstructed image');
end

if Bound==1 %Show the Fourier supports
   Show_Curvelets_boundaries(f,Bw,Bt,params.option);
end