function rec = iEWT2D_Mapping(ewtc,mfb)
%==========================================================================
% function rec = iEWT2D_Mapping(ewtc,mfb)
%
% This function performs the inverse Empirical Mapping Wavelet Transform,
% returning a reconstruction of the image.
%   
% Input:
%   -ewtc: empirical wavelet coefficients 
%   -mfb: Corresponding empirical wavelet filters
% Output:
%   -rec: reconstructed image
%
% Author: Charles Gerard-Lucas - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2024
% Version: 1.0
%==========================================================================
dual_sum = abs(ifftshift(mfb{1})).^2;
rec=fft2(ewtc{1}).*ifftshift(mfb{1}).*(abs(ifftshift(mfb{1}))>eps);
for i = 2:length(mfb)
    rec = rec + fft2(ewtc{i}).*ifftshift(mfb{i}).*(abs(ifftshift(mfb{i}))>eps);
    dual_sum = dual_sum + abs(ifftshift(mfb{i})).^2;
end

dual_sum(rec==0) = 1;
rec = ifft2(rec./dual_sum);