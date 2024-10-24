function Show_EWT2D_Filters(fil,shift)

%===========================================================
%
% function Show_EWT2D_Filters(fil,shift)
%
% This function permits to plot the Fourier magnitude of all 
% 2D EWT Littlewood-Paley filters.
%
% Input:
%   - fil: cell containing all the filters
%   - shift: 0=no fftshift ; 1=fftshift
%
% Author: J.Gilles
% Institution: UCLA - Department of Mathematics
% email: jegilles@math.ucla.edu
% Date: March, 1st, 2013
%
%===========================================================
p=ceil(length(fil)/2);

figure;

if shift==1
    for n=1:length(fil)
        subplot(2,p,n);imshow(fftshift(fil{n}),[]);
    end
else
    for n=1:length(fil)
        subplot(2,p,n);imshow(fil{n},[]);
    end

end