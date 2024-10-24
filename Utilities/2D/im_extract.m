function fe=im_extract(f,n)
%==========================================================================
% function fe = im_extract(f,n)
%
% This function removes edges of image f by n pixels.
%   
% Input:
%   -f: input image
%   -n: desired edge width 
% Output:
%   -fe: extracted image
%
% Author: Jerome Gilles
%==========================================================================

fe=f(n+1:end-n,n+1:end-n);