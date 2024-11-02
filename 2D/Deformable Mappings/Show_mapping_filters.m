function Show_mapping_filters(mfb)

%===========================================================
%
% function Show_mapping_filters(mfb)
%
% This function permits to plot the empirical wavelet filters.
%
% Input:
%   - mfb: structure of wavelet filters
%
% Author: C.-G. Lucas
% Institution: SDSU - Department of Mathematics
% Date: April, 19th, 2024
%
%===========================================================


p=ceil(length(mfb)/2);
figure; clf;
tiledlayout(2,p,'TileSpacing','tight','Padding','tight')
for k=1:length(mfb)
    nexttile
    imshow(abs(mfb{k})); axis tight; axis off; colormap gray; hold on;
    pbaspect([1 1 1])
end
