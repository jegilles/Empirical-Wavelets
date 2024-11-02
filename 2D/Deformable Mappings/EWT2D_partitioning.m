function [maxima,supports] = EWT2D_partitioning(f,params)

%==========================================================================
% function [maxima,supports] = EWT2D_partitioning(ff,params)
%
% Input:
%   -f: input image
%   -params: empirical wavelet parameters (see utilities)
%
% This function returns a partition of the Fourier domain following the
% Voronoi or Watershed method.
% Note: harmonic modes centered at the edge of the Fourier domain are discarded
% as they usually are not interested and/or of very small energy compared to the
% main modes. If you need to modofy this behavior, you can do so by changing the
% parameter params.edges whose value correspond to the width (in pixels) of the
% stripe around the image that is discarded.
%
% Output:
%   -maxima: detected local maxima representing modes
%   -supports: images containing the obtained partition
%
% Author: Charles-Gerard Lucas
% Institution: SDSU - Department of Mathematics
% Year: 2024
% Version: 1.0

ff = fftshift(fft2(f)); [h,w] = size(ff);
[absff,extH,extW] = Symmetrize_Fourier(abs(ff));
if params.log==1, absff=log(absff); end
absff = EWT_SpectrumRegularize_2D(absff,params);
maxima = EWWT_Get_Maxima(absff,params,extH,extW);

switch params.method
    case 'Watershed'        
        [supports,num_regions] = EWWT_Get_Boundaries(absff,maxima,params);
    case 'Voronoi'
        [~,vorcel] = EWT2D_Voronoi_Partition(maxima,size(absff));
        if isreal(f), vorcel = EWT2D_merge_symmetric(vorcel); end
        num_regions = 1:length(vorcel);
        supports=zeros(size(vorcel{1}));
        j = 0;
        for k=1:length(num_regions)
            j = j+1;
            if ~isempty(vorcel{k}), supports = supports+j*vorcel{k}; end
        end
end

% put null frequency in first position
approx = find(num_regions == supports(floor(h/2)+1,floor(w/2)+1));
supports_tmp = supports;
supports_tmp(supports==approx) = 1;
supports_tmp(supports==1) = approx;
supports = supports_tmp;

end