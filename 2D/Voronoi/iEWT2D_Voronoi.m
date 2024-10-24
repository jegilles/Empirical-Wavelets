function rec = iEWT2D_Voronoi(ewtc,mfb)
%==========================================================================
% function rec = iEWT2D_Voronoi(ewtc,mfb)
%
% This function performs the inverse Empirical Voronoi Wavelet Transform,
% returning a reconstruction of the image.
%   
% Input:
%   -ewtc: empirical wavelet coefficients 
%   -mfb: Corresponding empirical wavelet filters
%
% Output:
%   -rec: reconstructed image
%
% Author: Basile Hurat - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
%==========================================================================

dual_sum = mfb{1}.^2;
rec=fftshift(fft2(ewtc{1})).*mfb{1};
for i = 2:length(mfb)
    rec = rec + fftshift(fft2(ewtc{i})).*mfb{i};
    dual_sum = dual_sum + mfb{i}.^2;
end
rec = ifft2(ifftshift(rec./dual_sum));