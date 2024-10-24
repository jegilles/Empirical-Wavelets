function rec = iEWT2D_Watershed(ewtc,mfb)
%==========================================================================
% function rec = iEWT2D_Watershed(ewtc,mfb)
%
% This function performs the inverse Empirical Watershed Wavelet Transform,
% returning a reconstruction of the image.
%   
% Input:
%   -ewtc: empirical wavelet coefficients 
%   -mfb: Corresponding empirical wavelet filters
% Output:
%   -rec: reconstructed image
%
% Author: Basile Hurat - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
%==========================================================================
if nargin == 2
    imgReal = 1;
end
    
dual_sum = ifftshift(mfb{1}).^2;
rec=fft2(ewtc{1}).*ifftshift(mfb{1});
for i = 2:length(mfb)
    rec = rec + fft2(ewtc{i}).*ifftshift(mfb{i});
    dual_sum = dual_sum + ifftshift(mfb{i}).^2;
end

rec = ifft2(rec./dual_sum);