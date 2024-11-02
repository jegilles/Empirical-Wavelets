function Show_mapping_boundaries(f,supports,loga)
%==========================================================================
%function Show_mapping_boundaries(f,supports,loga)
%
% This function displays the detected boundaries (in red) and
% detected maxima in green) overtop the image's magnitude spectrum.
%
% Input:
%   -f: Image
%   -supports: Image containing the labelled partition
%   -loga: 1=logarithm of spectrum; 0 otherwise
%
% Author: Charles Lucas - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2024
% Version: 1.0
%==========================================================================

% remove the extensions
extH=0; extW=0;
if rem(size(f,1),2)==0, extH=1; end; if rem(size(f,2),2)==0, extW=1; end
if extH==1, supports = supports(1:end-1,:); end
if extW==1, supports = supports(:,1:end-1); end

if loga == 1
    ff = log(1+abs(fftshift(fft2(f))));
else 
    ff = abs(fftshift(fft2(f)));   
end
ff=ff/max(ff(:));

figure;
imshow(labeloverlay(ff,boundarymask(supports),'Transparency',0,'Colormap','hot'))