function [centers, plane] = EWWT_Get_Maxima(absff,params,extH,extW)
%==========================================================================
% function centers = EWWT_Get_Maxima(absff,params,extH,extW)
%
% This function uses the EWT_GSS_2D function to build a scale-space 
% representation of f via successive separable convolutions with a 
% discrete Gaussian Kernel. Then, it post-processes the detected maxima.
% Note: harmonic modes centered at the edge of the Fourier domain are discarded
% as they usually are not interested and/or of very small energy compared to the
% main modes. If you need to modofy this behavior, you can do so by changing the
% parameter params.edges whose value correspond to the width (in pixels) of the
% stripe around the image that is discarded.
%
% Inputs:
% -absff: the input image
% -params: empirical wavelet parameters (see utilities folder)
% -extH: flag to say wether original image's HEIGHT was even and had to
%   be extended to get symmetric spectrum
% -extW: flag to say wether original image's WIDTH was even and had to
%   be extended to get symmetric spectrum
%
% Outputs:
% -centers: collection of the coordinates of the significant maxima
%
% Author: Basile Hurat - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
% =========================================================================
[hh,ww] = size(absff);
[centers,~,plane] = EWT_GSS_2D(absff,params);

% Discard centers detected with params.edges pixels from the edge of the Fourier
% domain.
% If edges parameter greater than zero, set that much of the edge equal to 0
if params.edges > 0 || params.includeCenter == 1
    center_flag = 1; %Assume we don't have center
    i = 1;
    while i <= size(centers,1)
        if centers(i,1) == round(hh/2) && centers(i,2) == round(ww/2)
            center_flag = 0; %Change flag if we find center
        end
        if centers(i,1) <= params.edges || centers(i,2) <= params.edges
            centers(i,:) = [];
        elseif centers(i,1) >= hh+1-params.edges || centers(i,2) >= ww+1-params.edges
                centers(i,:) = [];        
        else
            i = i + 1;
        end        
    end
    if center_flag == 1
        centers =[centers; [round(hh/2),round(ww/2)]];
    end
end

%If image was extended due to being even, remove those extensions
if extH == 1
    i = 1; 
    while i <= size(centers,1)
        if centers(i,1) == 1 || centers(i,1) == hh
            centers(i,:) = [];
        end
        i = i + 1;
    end
end
if extW == 1
    i = 1; 
    while i <= size(centers,1)
        if centers(i,2) == 1 || centers(i,2) == ww
            centers(i,:) = [];
        else
            i = i + 1;
        end
    end
end