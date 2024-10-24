function fe = im_fourier_extend(f,n)
%==========================================================================
% function fe = im_fourier_extend(f,n)
%
% This function extends the dimensions of image f by n pixels on all sides. 
% It periodize the image.
%   
% Input:
%   -f: input image
%   -n: desired pixel extension
% Output:
%   -fe: extended image
%
% Author: Jerome Gilles
%==========================================================================

fe=zeros(size(f,1)+2*n,size(f,2)+2*n);

% copy f in the center of the extension
fe(n+1:n+size(f,1),n+1:n+size(f,2))=f;

% extend horizontally
fe(1:n,:)=fe(end-2*n+1:end-n,:);
fe(end-n+1:end,:)=fe(n+1:2*n,:);
 

% extend vertically
fe(:,1:n)=fe(:,end-2*n+1:end-n);
fe(:,end-n+1:end)=fe(:,n+1:2*n);
