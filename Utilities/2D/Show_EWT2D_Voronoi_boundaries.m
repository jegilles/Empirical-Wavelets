function Show_EWT2D_Voronoi_boundaries(f,vor,color,logspec)

%===============================================================================
% function Show_EWT2D_Voronoi_boundaries(f,vor,color,logspec)
%
% This function plot the edges of the Voronoi partition onto the magnitude of
% the Fourier spectrum of the input image.
% 
% Inputs:
%   -f: input image
%   -vor: image containong the plot of the Voronoi partition
%   -color: vector corresponding to the wanted color for the partition edges.
%   This must be a vector of the form [R G B] where the quantities R,G,B are
%   provided within [0,1]. If this parameter is not provided, the plot will be
%   made in red by default ([1 0 0]).
%   -logspec: flag to indicate if the function plots the logarithm of the image
%   spectrum (logspec = 1). If not provided, the default behavior is to NOT plot
%   the logarithm (i.e logspec = 0).
%
% Author: Jerome Gilles
% Institution: Dept of Mathematics & Statistics, San Diego State University
% Version: 1.0
% Year: 2022
%===============================================================================

if nargin < 3 || isempty(color)
    color = [1 0 0];
end

if nargin < 4
    logspec = 0;
end

% Find the edges of the Voronoi partition
gr=zeros(size(f));
gr(1:end-1,:) = gr(1:end-1,:)+diff(vor,1,1);
gr(:,1:end-1) = abs(gr(:,1:end-1))+abs(diff(vor,1,2));

% Get the magnitude of the Fourier spectrum of f
absff = abs(fftshift(fft2(f)));
if logspec == 1
    absff = log(1+absff);
end


%normalize between [0,1]
absff=(absff-min(absff(:)))/(max(absff(:))-min(absff(:)));
% make it an RGB image
f=repmat(absff, [1 1 3]);

% tag each edge pixel to the wanted color
for i=1:size(gr,1)
    for j=1:size(gr,2)
        if gr(i,j)~=0
            f(i,j,:) = color;
        end
    end
end

%imwrite(f,'vorospecspec.png');

figure;
imshow(f,[]);

