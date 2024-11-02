%% This script allows you to test the empirical wavelet transform using mappings
% It shows the usage of the functions as well as demonstrating some of the
% results found in C.-G. Lucas and J. Gilles, "Multidimensional Empirical
% Wavelet Transform" (2024)
clear all
close all

%% Signal and results to display
% Choose an image from list (Texture, Lena, Barb) or load in your own to
% variable f
signal = "texture";

% Choose the functions to use
Partitioning = 1;   % Run partitioning outside the main function
Mapping = 1;  % Run mapping estimation outside the main function

switch lower(signal)
    case 'texture'
        load('texture.mat');
    case 'lena'
        load('lena.mat');
    case 'barb'
        load('barb.mat');
end

f=f(1:128,1:128);

%% Partitioning parameter setup
% Begin by running default parameters in case you miss a parameter
params = EWTDefaultParams(); 

% Choose a partitioning method
params.method = 'Watershed'; % 'Voronoi' or 'Watershed', or define params.supports manually

% Number of edges to remove
params.edges = 2;

% Choose the wanted regularization (none,gaussian,average,closing)
params.reg = 'none';
params.lengthFilter = 10;
params.sigmaFilter = 1.5;

% Choose the wanted detection method for the scales (scalespace or nlevel)
params.detect = 'scalespace'; 
params.typeDetect='otsu'; %for scalespace:otsu,halfnormal,empiricallaw,kmeans

% Perform the detection on the log spectrum instead the spectrum
params.log = 1;

% Scale-space step-size
params.t = .8;

%% Wavelet filterbank parameter setup
% Choose a mother wavelet
name = 'Disk'; % 'Gabor', 'Shannon', 'Disk' or 'Square'
params_wavelet.tau = .2; % for mother wavelet 'Disk' and 'Square'
[wavelet.Lambda,wavelet.mother_wavelet] = EWT2D_mother_wavelet(f,name,params_wavelet);
wavelet.normalization = 'none'; % 'none' or 'det'

% Choose demons algorithm
algo = 'diffeomorphic';
[h,w] = size(f); max_level = floor(log2(min(h,w)));
switch algo
    case 'Thirion'     
        demons.algo = 'Thirion';
        demons.level = max_level-1:max_level;
        demons.diffus = .3:.01:.5;
    case 'diffeomorphic'
        demons.algo = 'Vercauteren';
        demons.field = 'exp';
        demons.level = max_level-1:max_level;
        demons.diffus = .3:.01:.5;
        demons.spatial_uncertainty = 5;
    case 'additive'
        demons.algo = 'Vercauteren';
        demons.field = 'identity';
        demons.level = max_level-1:max_level;
        demons.diffus = .3:.01:.5;
        demons.spatial_uncertainty = 4;
end
mapping.demons = demons;

%% Partitioning
if Partitioning==1
   [partition.maxima,partition.supports] = EWT2D_partitioning(f,params);
   Show_mapping_boundaries(f,partition.supports,params.log);
else
   partition = params; % run in EWT2D_Mapping
end

%% Mapping estimation
if Mapping==1
    mapping = EWT2D_demons_mappings(wavelet.Lambda,partition.supports,demons);
    Show_mappings(partition.supports,mapping.phi_Omega,~isreal(f))
else
    mapping.demons = demons; % run in EWT2D_Mapping
end

%% Empirical wavelet transform
[ewtc, mfb, centers, supports, phi_Omega, phi] = EWT2D_Mapping(f,partition,wavelet,mapping);

%% Show filterbank and wavelet spectrum
Show_mapping_filters(mfb)
Show_EWT2D_mapping_spectrum(ewtc);

%% Reconstruction
figure;subplot(1,2,1);imshow(f,[]);title('Original image');
rec = iEWT2D_Mapping(ewtc,mfb);
fprintf('Reconstruction error: %e\n',norm(f(:)-rec(:),'fro'));
subplot(1,2,2);imshow(real(rec),[]);title('Reconstructed image');
