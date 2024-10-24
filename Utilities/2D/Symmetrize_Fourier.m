function [sym,extH,extW] = Symmetrize_Fourier(in)

%===============================================================================
% function sym = Symmetrize_Fourier(in)
%
% This function returns a symetrized Fourier transform by properly extending
% the expected directions. It returns a new spectrum where each dimension is
% odd
%
% Inputs:
%   -in: magnitude of the spectrum to symetrize
%
% Outputs:
%   -sym: symetrized spectum
%   -extH: tag: 1 if there is vertical extension, 0 otherwise
%   -extW: tag: 1 if there is horizontal extension, 0 otherwise
%
%
% Author: Jerome Gilles
% Institution: San Diego State University - Dept of Mathematics & Statistics
% Version: 1.0
% Year: 2022
%===============================================================================

if (rem(size(in,1),2)==0) && (rem(size(in,2),2)~=0) % vertical extension
    sym = zeros(size(in,1)+1,size(in,2));
    sym(1:end-1,:) = in;
    sym(end,:) = in(1,end:-1:1);
    extH=1 ; extW=0;
elseif (rem(size(in,1),2)~=0) && (rem(size(in,2),2)==0) % horizontal extension
    sym = zeros(size(in,1),size(in,2)+1);
    sym(:,1:end-1) = in;
    sym(:,end) = in(end:-1:1,1);
    extH=0 ; extW=1;
elseif (rem(size(in,1),2)==0) && (rem(size(in,2),2)==0) %extension in both directions
    sym = zeros(size(in,1)+1,size(in,2)+1);
    sym(1:end-1,1:end-1) = in;
    sym(end,2:end-1) = in(1,end:-1:2);
    sym(2:end-1,end) = in(end:-1:2,1);
    sym(1,1) = in(1,1);
    sym(1,end) = in(1,1);
    sym(end,1) = in(1,1);
    sym(end,end) = in(1,1);
    extH=1 ; extW=1;
else % already symetric
    sym=in;
    extH=0 ; extW=0;
end
