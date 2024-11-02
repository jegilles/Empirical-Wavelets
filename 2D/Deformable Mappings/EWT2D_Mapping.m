function [ewtc,mfb,centers,supports,phi_Omega,phi] = EWT2D_Mapping(f,params,wavelet,mapping)

%==========================================================================
% function [ewtc,mfb,centers,supports,phi_Omega,phi] = EWT2D_Mapping(f,params)
%
% This function performs the 2D Empirical Wavelet Transform.
% The partitioning is achieved by a combination of scale-space
% representations for mode detection and the Watershed or Voronoi
% partitioning. Then, the empirical wavelet system is built by estimating
% mappings from the selected mother wavelet's Fourier support to the 
% partition supports. The whole procedure is detailed in C.-G. Lucas and 
% J. Gilles, "Multidimensional Empirical Wavelet Transform" (2024)
%   
% Input:
%   -f: input image
%   -params: empirical wavelet parameters (see utilities)
%
% Output:
%   -ewtc: cell containing each filtered output subband (ewtc{1} is the 
%   lowpass subband which corresponds to the region that contains the 
%   origin and ewtc{k} are the bandpass filtered images corresponding to 
%   the other regions
%   -mfb: cell containing the set of empirical filters in the Fourier 
%   domain (the indexation is the same as ewtc above)
%   -centers: detected local maxima representing modes
%   -supports: images containing the obtained partition
%   -phi_Omega: images containing the estimated partition
%   -phi: displacement field which corresponds to the estimated mapping
%
% Author: Charles-Gerard Lucas - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2024
% Version: 1.0
%==========================================================================

params.complex = ~isreal(f);
[h,w] = size(f); if rem(h,2)==0, extH=1; end; if rem(w,2)==0, extW=1; end

% Get Fourier transform and its partition
ff = fftshift(fft2(f));
if isfield(params,'supports')
    supports = params.supports; maxima = params.maxima;
else
    [maxima,supports] = EWT2D_partitioning(f,params);
end

% Create the filterbank
[mfb,phi_Omega,phi] = EWT2D_Mapping_Filterbank(supports,wavelet,mapping,extH,extW,params.complex);

% remove the extensions
if extH==1, supports = supports(1:end-1,:); end
if extW==1, supports = supports(:,1:end-1); end

% get centers
centers = cell([round(length(maxima)/2),1]);
for i = 1:length(centers)
    if maxima(i,1) == floor(h/2)+1 && maxima(i,2) == floor(w/2)+1
        centers{i} = maxima(i,:);
    else
        centers{i} = [maxima(i,:);maxima(end-i+1,:)];
    end
end

% Filter the input image and return EWT coefficients
ewtc = cell(size(mfb));
for k = 1:length(ewtc)
    ewtc{k} = ifft2(ifftshift(ff.*conj(mfb{k})));
end