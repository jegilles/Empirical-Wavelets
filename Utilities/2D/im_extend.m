function fe = im_extend(f,n)
%==========================================================================
% function fe = im_extend(f,n)
%
% This function extends the dimensions of image f by n pixels on all sides. 
% It mirrors the edges.
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
fe(n+1:n+size(f,1),n+1:n+size(f,2))=f;
fe(1:n,:)=fe(2*n:-1:n+1,:);
fe(end:-1:end-n+1,:)=fe(end-2*n:end-n-1,:);
fe(:,1:n)=fe(:,2*n:-1:n+1);
fe(:,end:-1:end-n+1)=fe(:,end-2*n:end-n-1);