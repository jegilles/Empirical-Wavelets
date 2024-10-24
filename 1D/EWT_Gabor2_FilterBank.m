function [mfb,ws]=EWT_Gabor2_FilterBank(boundaries,N,cpx)

% =========================================================================
% function mfb=EWT_Gabor2_FilterBank(boundaries,N,cpx)
%
% This function generate the Gabor filter bank (scaling function
% + wavelets) corresponding to the provided set of frequency segments.
% In this version, the highpass filters are made of half a Gaussian and
% half a constant.
% For the real case, the filters are defined on symmetric supports while
% on each individual support for the complex case.
%
% Input parameters:
%   -boundaries: vector containing the boundaries of frequency segments (0
%                and pi must NOT be in this vector)
%   -N: signal length
%   -cpx: 1 if the original signal was complex, 0 otherwise
%
% Output:
%   -mfb: cell containing each filter (in the Fourier domain), the scaling
%         function comes first and then the successive wavelets (for the
%         complex case, filters in the positive frequencies come before
%         filters in the negative frequencies)
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics & Statistics
% Year: 2020
% Version: 1.0
% =========================================================================

Npic=length(boundaries);

%create the frequency axis [0,2pi)
Mi=floor(N/2)+1;
ws=(0:2*pi/N:2*pi-2*pi/N)'; 
ws(Mi+1:end)=-2*pi+ws(Mi+1:end); %switch second half to negative frequencies


if cpx == 0
    mfb=cell(Npic+1,1);

    % scaling function
    mfb{1}=zeros(length(ws),1);
    om=2*boundaries(1);
    for k=1:length(mfb{1})
      mfb{1}(k)=exp(-pi*(2.5*abs(ws(k))/om)^2)/sqrt(om);
    end
    
    % each wavelet except the last one
    for l=1:Npic-1
       mfb{l+1}=zeros(length(ws),1);

        om=boundaries(l+1)-boundaries(l);
        wm=(boundaries(l+1)+boundaries(l))/2;
        for k=1:length(mfb{l+1})
          mfb{l+1}(k)=exp(-pi*(2.5*(abs(ws(k))-wm)/om)^2)/sqrt(om);
        end
    end
    
    % last one
    mfb{end}=zeros(length(ws),1);

    om=pi-boundaries(end);
    wm=(boundaries(end)+pi)/2;

    for k=1:length(mfb{end})
        if abs(ws(k))<wm
            mfb{end}(k)=exp(-pi*(2.5*(abs(ws(k))-wm)/om)^2)/sqrt(om);
        else
            mfb{end}(k)=1/sqrt(om);
        end
    end
else
    boundaries(boundaries>=pi)=-2*pi+boundaries(boundaries>=pi);
    % note that, contrarily to the papers, the negative boundaries are
    % listed AFTER the positive ones.
    
    mfb=cell(Npic,1);

    % scaling function
    mfb{1}=zeros(length(ws),1);

    om=boundaries(1)-boundaries(end);
    wm=(boundaries(1)+boundaries(end))/2;
    for k=1:length(mfb{1})
      mfb{1}(k)=exp(-pi*(2.5*(ws(k)-wm)/om)^2)/sqrt(om);
    end
    
    % each wavelet (High frequencies in a separate way)
    shift=1;
    for l=1:Npic-1
        if boundaries(l)>boundaries(l+1) %high pass filters
            %positive frequencies
            mfb{l+shift}=zeros(length(ws),1);
            om=pi-boundaries(l);
            wm=(boundaries(l)+pi)/2;
            for k=1:length(mfb{1})
                if ws(k)<wm
                    mfb{l+shift}(k)=exp(-pi*(2.5*(ws(k)-wm)/om)^2)/sqrt(om);
                else
                    mfb{l+shift}(k)=1/sqrt(om);
                end
            end
            shift=2;
            %negative frequencies
            mfb{l+shift}=zeros(length(ws),1);
            om=abs(pi+boundaries(l+1));
            wm=(-pi+boundaries(l+1))/2;
            for k=1:length(mfb{1})
                if wm<=ws(k)
                    mfb{l+shift}(k)=exp(-pi*(2.5*(ws(k)-wm)/om)^2)/sqrt(om);
                else
                    mfb{l+shift}(k)=1/sqrt(om);
                end
            end
        else %individual wavelets
            mfb{l+shift}=zeros(length(ws),1);
            
            om=boundaries(l+1)-boundaries(l);
            wm=(boundaries(l)+boundaries(l+1))/2;
            for k=1:length(mfb{1})
              mfb{l+shift}(k)=exp(-pi*(2.5*(ws(k)-wm)/om)^2)/sqrt(om);
            end
        end
    end
end