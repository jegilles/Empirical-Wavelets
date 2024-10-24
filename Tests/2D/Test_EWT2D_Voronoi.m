%% This script allows you to test the empirical Voronoi wavelet transform

%% Parameter setup
% Choose an image from list (Texture, Lena, Barb) or load in your own to
% variable f
signal = "Texture";

% Choose the results you want to display (Show=1, Not Show=0)
Bound = 1;   % Display the detected boundaries on the spectrum
Comp = 0;    % Display the EWT components
Rec = 0;     % Display the reconstructed signal

%Begin by running default parameters in case you miss a parameter
params = EWTDefaultParams(); 

% Choose the wanted regularization (none,gaussian,average,closing)
params.reg = 'none';
params.lengthFilter = 10;
params.sigmaFilter = 1.5;


% Choose the wanted detection method for the scales (scalespace or nlevel)
params.detect = 'scalespace'; 
params.typeDetect='otsu'; %for scalespace:otsu,halfnormal,empiricallaw,mean,kmeans

% Perform the detection on the log spectrum instead the spectrum
params.log=0;

switch lower(signal)
    case 'texture'
        load('texture.mat');
        %f=noisy;
    case 'lena'
        load('lena.mat');
    case 'barb'
        load('barb.mat');
end

%% We perform the empirical watershed wavelet transform on the signal f
[ewtc, mfb, centers, vorpartition,plane] = EWT2D_Voronoi(f,params);

%% We show desired results
if Comp==1 %Show watershed wavelet components
   Show_EWT2D(ewtc,params);
end

if Rec==1 %Show the reconstructed image
   figure;subplot(1,2,1);imshow(f,[]);title('Original image');
   rec = iEWT2D_Voronoi(ewtc,mfb);
   fprintf('Reconstruction error: %e\n',norm(f(:)-rec(:),2));
   subplot(1,2,2);imshow(rec,[]);title('Reconstructed image');
end

if Bound==1 %Show the Fourier supports
   Show_EWT2D_Voronoi_boundaries(f,vorpartition,[],1)
end
