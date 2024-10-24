function [ewtc,mfb,centers,supports] = EWT2D_Watershed(f,params)
%==========================================================================
% function [ewtc,mfb,centers,supports] = EWT2D_Watershed(f,params)
%
% This function performs the Empirical Watershed Wavelet Transform.
% The partitioning is achieved by a combination of scale-space
% representations for mode detection and the watershed transform for actual
% partitioning. The watershed transform is is performed by the built-in
% matlab function based on F.Meyer's watershed algorithm. After, it uses
% the arbitrary shape empirical wavelet transform. 
% You can read the details in Hurat, Alvarado, and Gilles "The Empirical 
% Watershed Wavelet Transform" (J. Imaging 2020)
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
%
% Author: Basile Hurat - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
%==========================================================================


%Get Fourier transform and its magnitude
[h,w] = size(f);

if isreal(f)
    params.complex = 0;
else
    params.complex = 1;
end

ff = fftshift(fft2(f));
ff_orig = ff;
extH = 0; extW = 0;

%Make the image odd (to guarantee symmetry)
if mod(h,2) == 0
    ff(h+1,:) = ff(1,:);
    extH = 1;
end
if mod(w,2) == 0
    ff(:,w+1) = ff(:,1);
    extW = 1;
end

absff = abs(ff);
if params.log==1
    absff=log(absff);
end

% Global trend removal and regularization
presig = absff; %EWT_RemoveTrend_2D(absff,params); no trend removal at this time
absff = EWT_SpectrumRegularize_2D(presig,params);


%Get meaningful maxima, make cell for centers
maxima = EWWT_Get_Maxima(absff,params,extH,extW);

centers = cell([round(length(maxima)/2),1]);
for i = 1:length(centers)
    if maxima(i,1) == floor(h/2)+1 && maxima(i,2) == floor(w/2)+1
        centers{i} = maxima(i,:);
    else
        centers{i} = [maxima(i,:);maxima(end-i+1,:)];
    end
end


% Use watershed to partition the Fourier domain starting with the detected 
% maxima as flooding points
[supports,num_regions] = EWWT_Get_Boundaries(absff,maxima,params);


% Create the filter bank
mfb = EWT2D_Watershed_Filterbank(supports,num_regions,extH,extW,params.tau);

% Switch the filter including the zero frequency as the first showing up in the
% cell
approx = find(num_regions == supports(floor(h/2)+1,floor(w/2)+1));
temp = mfb{1};
mfb{1} = mfb{approx};
mfb{approx} = temp;
temp = centers{1};
centers{1} = centers{approx};
centers{approx} = temp;


% remove the extensions
if extH==1
    supports = supports(1:end-1,:);
end
if extW==1
    supports = supports(:,1:end-1);
end

% Filter the input image and return EWT coefficients
ewtc = cell([length(mfb) 1]);
if params.complex == 0
    for k = 1:length(mfb)
        ewtc{k} = real(ifft2(ifftshift(ff_orig.*conj(mfb{k}))));
    end
else
    for k = 1:length(mfb)
        ewtc{k} = ifft2(ifftshift(ff_orig.*conj(mfb{k})));
    end
end