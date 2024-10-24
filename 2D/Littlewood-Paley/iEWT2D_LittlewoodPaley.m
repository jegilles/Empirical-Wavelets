function rec=iEWT2D_LittlewoodPaley(ewt,mfb)

% ======================================================================
% function rec=iEWT2D_LittlewoodPaley(ewt,mfb)
% 
% Perform the Inverse 2D Littlewood-Paley Empirical Wavelet Transform of 
% ewt accordingly to the filter bank mfb
%
% Inputs:
%   -ewt: cell containing the 2D EWT components
%   -mfb: filter bank used during the EWT
%
% Output:
%   -rec: reconstructed image
%
% Author: Jerome Gilles - Giang Tran
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
% ======================================================================

%We perform the adjoint operator to get the reconstruction
for k=1:length(ewt)
    if k==1
        rec=zeros(size(ewt{1}));
    end
    rec=rec+real(ifft2(fft2(ewt{k}).*mfb{k}));
end