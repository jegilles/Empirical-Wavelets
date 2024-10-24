function EWT_Plot_Shannon_FilterBank(mfb)

%==========================================================================
% function EWT_Plot_Shannon_FilterBank(mfb,cpx)
%
% This function plots the Meyer Complex filter bank (both mgnitude and
% phase
%
% Input:
%   -mfb: Meyer complex filter bank
%
%==========================================================================

Mi=floor(length(mfb{1})/2)+1;
xi=(0:2*pi/length(mfb{1}):2*pi-2*pi/length(mfb{1}))'; 
xi(Mi+1:end)=-2*pi+xi(Mi+1:end); %switch second half to negative frequencies
xi=fftshift(xi);

mag=cell(size(mfb));
ang=cell(size(mfb));

for k=1:length(mfb)
    mag{k}=fftshift(abs(mfb{k}));
end

for k=1:length(mfb)
     ang{k}=fftshift(angle(mfb{k}));
end

for k=1:length(mfb)
    subplot(2,1,1);plot(xi,mag{k}); title('Magnitude of wavelet filters');
    if k==1
        hold on
    end
    subplot(2,1,2);plot(xi,ang{k}); title('Phase of wavelet filters');
    if k==1
        hold on
    end
end
hold off