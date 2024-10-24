function [ewtc,mfb,maxima,vorpartition,plane] = EWT2D_Voronoi(f,params)

%==========================================================================
% function [ewtc,mfb,maxima,vorpartition] = EWT2D_Voronoi(f,params)
%
% This function compute the 2D EWT based on a Voronoi partitioning of the
% Fourier domain. If the real transform is required (i.e params.complex ~= 1),
% symmetric cells are grouped together to guarantee symmetric filters.
% Note: harmonic modes centered at the edge of the Fourier domain are discarded
% as they usually are not interested and/or of very small energy compared to the
% main modes. If you need to modofy this behavior, you can do so by changing the
% parameter params.edges whose value correspond to the width (in pixels) of the
% stripe around the image that is discarded.
%   
% Input:
%   -f: input image
%   -params: transform parameters
%
% Output:
%   -ewtc: collection of outputs of each EW filter
%   -mfb: the built filter bank
%   -maxima: coordinates of each meaningful detected maxima
%   -vorpartition: detected Voronoi partition
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
%==========================================================================

%Get Fourier transform and its magnitude
ff_orig = fftshift(fft2(f));

% symmetrize the spectrum if needded
[absff,extH,extW] = Symmetrize_Fourier(abs(ff_orig));

if params.log==1
    absff=log(absff);
end

% regularization (if requested)
absff = EWT_SpectrumRegularize_2D(absff,params);

%Get meaningful maxima, make cell for centers
[maxima,plane] = EWWT_Get_Maxima(absff,params,extH,extW);

% %build the voronoi cells
[vorpartition,vorcel] = EWT2D_Voronoi_Partition(maxima,size(absff));
vorpartition = UnSymmetrize_Fourier(vorpartition,extH,extW);

if params.complex ~= 1 %group symetric cells to guarantee a real transform
    vorocells = EWT2D_merge_symmetric(vorcel);
    % Create the filterbank
    mfb = EWT2D_Voronoi_Filterbank(size(absff),vorocells,params.tau,extH,extW);
else
    mfb = EWT2D_Voronoi_Filterbank(size(absff),vorcel,params.tau,extH,extW);
end

% Perform the filtering
ewtc = cell(size(mfb));

for i=1:length(ewtc)
    ewtc{i}=ifft2(ifftshift(ff_orig.*mfb{i})); 
end