function Show_Watershed_boundaries(f,supports,loga)
%==========================================================================
%function Show_Watershed_boundaries(f,supports,loga)
%
% This function displays the detected watershed boundaries (in red) and
% detected maxima in green) overtop the image's magnitude spectrum.
%
% Input:
%   -f: Image
%   -supports: Image containing the labelled partition
%   -loga: 1=logarithm of spectrum; 0 otherwise
%
% Author: Basile Hurat - Jerome Gilles
% Institution: SDSU - Department of Mathematics
% Year: 2022
% Version: 1.0
%==========================================================================
if loga == 1
    ff = log(1+abs(fftshift(fft2(f))));
else 
    ff = abs(fftshift(fft2(f)));    
end
ff=ff/max(ff(:));

figure;
imshow(labeloverlay(ff,boundarymask(supports),'Transparency',0,'Colormap','hot'))