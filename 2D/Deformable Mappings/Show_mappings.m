function Show_mappings(supports,phi_Omega,cpx)

%===========================================================
%
% function Show_mapping(ewt)
%
% This function permits to plot the mapping estimates of the Fourier
% partition of the EWT2D
%
% Input:
%   - supports: matrix of the Fourier supports
%   - phi_Omega: preimage of the mother wavelet's support by the estimated
%   diffeomorphisms
%
% Author: C.-G. Lucas
% Institution: SDSU - Department of Mathematics
% Date: June, 28th, 2024
%
%===========================================================
p=ceil(length(phi_Omega)/2);
[h,w] = size(supports); Mi=floor(h/2)+1; Mj=floor(w/2)+1;

linewidth = 1.5;
figure; clf;
tiledlayout(2,p,'TileSpacing','tight','Padding','tight')
for k=1:length(phi_Omega)
    nexttile
    Omegak = (bwlabel((supports==k),8)==1);
    boundary = bwboundaries(Omegak,'noholes');
    switch cpx
        case 1
            imshow(phi_Omega{k}); axis tight; axis off; colormap gray; hold on;
            plot(boundary{1}(:,2), boundary{1}(:,1), 'r', 'LineWidth', linewidth);
        case 0
            if Omegak(Mi,Mj)
                imagesc(phi_Omega{k}); axis tight; axis off; colormap gray; hold on;
                plot(boundary{1}(:,2), boundary{1}(:,1), 'r', 'LineWidth', linewidth);
            else
                imagesc(phi_Omega{k}+flip(flip(phi_Omega{k},1),2)); axis tight; axis off; colormap gray; hold on;
                plot(boundary{1}(:,2), boundary{1}(:,1), 'r', 'LineWidth', linewidth);
                boundary2 = bwboundaries(flip(flip(Omegak,1),2),'noholes');
                plot(boundary2{1}(:,2), boundary2{1}(:,1), 'r', 'LineWidth', linewidth);
            end
    end        
    pbaspect([1 1 1])
end