function sym = EWT2D_Symmterize_Curvelet_Filter(fil)

%==========================================================================
% function sym = EWT2D_Symmterize_Curvelet_Filter(fil)
%
% This function returns a mirrored version of the input.
%
% Input:
%   fil: input image. IT MUST BE OF ODD SIZE!
%
% Output:
%   sym: mirrored version of the input
%
% Author: Jerome Gilles
% Institution: SDSU Dept of Mathematics & Statistics
% Version: 1.0 (2019)
%==========================================================================

sym=fil(end:-1:1,end:-1:1);