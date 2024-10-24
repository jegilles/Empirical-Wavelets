function Show_EWT_Filters(mfb,f)

% ====================================================================
% function Show_EWT_Filters(mfb,f)
%
% This function plots each filter on top of the magnitude spectrum of f.
% The plots are on the frequency interval [-pi,pi).
% NOTE: the magnitude spectrum has been normalized!
%
% Inputs:
%   -mfb: EWT wavelet filters
%   -f: input signal
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics & Statistics
% Year: 2019
% Version: 1.0
% =====================================================================

figure;
x=0:1/length(mfb{1}):(length(mfb{1})-1)/length(mfb{1});
x=2*pi*x-pi;

%don't plot more than 6 filters per figure
l=1;
if length(mfb)>6
    lm=6;
else
    lm=length(mfb);
end

ff=fftshift(abs(fft(f)));
ff=ff/max(ff);

for k=1:length(mfb)
    hold on; 
    if isreal(mfb{k})
        subplot(lm,1,l); plot(x,ff,x,fftshift(mfb{k}));
        title(['mfb(',num2str(k),')'])
        if mod(k,6) == 0
            figure;
            l=1;
        else
            l=l+1;
        end
    else
        subplot(lm,1,l); plot(x,ff,x,fftshift(real(mfb{k})));
        title(['real part of mfb(',num2str(k),')'])
        subplot(lm,2,l+1); plot(x,ff,x,fftshift(imag(mfb{k})));
        title(['imaginary part of mfb(',num2str(k),')'])
        if mod(k,6) == 0
            figure;
            l=1;
        else
            l=l+2;
        end
    end
end