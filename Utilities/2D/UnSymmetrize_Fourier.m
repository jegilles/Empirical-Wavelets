function usym = UnSymmetrize_Fourier(in,extH,extW)

%===============================================================================
% function usym = UnSymmetrize_Fourier(in,extH,extW)
%
% This function removes the extension of the spectrum that were introduced to
% symmetrize the spectrum.
%
% Inputs:
%   -in: magnitude of the spectrum to symetrize
%   -extH: tag: 1 if there is vertical extension, 0 otherwise
%   -extW: tag: 1 if there is horizontal extension, 0 otherwise
%
% Outputs:
%   -sym: original size spectum
%
%
% Author: Jerome Gilles
% Institution: San Diego State University - Dept of Mathematics & Statistics
% Version: 1.0
% Year: 2022
%===============================================================================

if (extH==1) && (extW==0) % remove vertical extension
    usym = in(1:end-1,:);
elseif (extH==0) && (extW==1) % remove horizontal extension
    usym = in(:,1:end-1);
elseif (extH==1) && (extW==1) %remove in both directions
    usym = in(1:end-1,1:end-1);
else % no extension to remove
    usym=in;
end
