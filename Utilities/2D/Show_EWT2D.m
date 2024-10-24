function Show_EWT2D(ewt,params)

%===========================================================
%
% function Show_EWT2D(ewt)
%
% This function permits to plot the outputs of the 2D EWT.
%
% Input:
%   - ewt: cell containing the EWT outputs
%
% Author: J.Gilles
% Institution: UCLA - Department of Mathematics
% email: jegilles@math.ucla.edu
% Date: March, 1st, 2013
%
%===========================================================
p=ceil(length(ewt)/2);


if params.complex ~= 1
    figure;
    for n=1:length(ewt)
        subplot(2,p,n);imshow(ewt{n},[]);
        name = sprintf('coef%d.png',n);
        imwrite((ewt{n}-min(ewt{n}(:)))/(max(ewt{n}(:))-min(ewt{n}(:))),name);
    end
else
    figure;
    for n=1:length(ewt)
        subplot(2,p,n);imshow(real(ewt{n}),[]);
    end
    sgtitle('Real parts');
    figure;
    for n=1:length(ewt)
        subplot(2,p,n);imshow(imag(ewt{n}),[]);
    end
    sgtitle('Imaginary parts');
end