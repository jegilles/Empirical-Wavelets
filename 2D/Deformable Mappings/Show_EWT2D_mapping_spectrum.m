function Show_EWT2D_mapping_spectrum(ewtc)

%===========================================================
%
% function Show_EWT2D_mapping_spectrum(ewtc)
%
% This function permits to plot the outputs of the 2D EWT.
%
% Input:
%   - ewt: cell containing the EWT outputs
%
% Author: C.-G. Lucas
% Institution: SDSU - Department of Mathematics
% Date: April, 19th, 2024
%
%===========================================================
p=ceil(length(ewtc)/2);

figure; clf;
tiledlayout(2,p,'TileSpacing','tight','Padding','tight')
for n=1:length(ewtc)
    nexttile; 
    imagesc(abs(ewtc{n}).^2); axis tight; axis off; colormap gray
    pbaspect([1 1 1])
end
