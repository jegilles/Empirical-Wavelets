function [mfb,ws]=EWT_Shannon_FilterBank(boundaries,N,cpx)

% =========================================================================
% function mfb=EWT_Shannon_FilterBank(boundaries,N,cpx)
%
% This function generate the Shannon filter bank (scaling function
% + wavelets) corresponding to the provided set of frequency segments.
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

    for k=1:length(mfb{1})
       if abs(ws(k))<boundaries(1)
          mfb{1}(k)=1;
       end
    end
    
    % each wavelet except the last one
    for l=1:Npic-1
       mfb{l+1}=zeros(length(ws),1);

        wl=boundaries(l+1)-boundaries(l);
        for k=1:length(mfb{l+1})
           if (boundaries(l)<=abs(ws(k))) && (abs(ws(k))<boundaries(l+1))
              mfb{l+1}(k)=exp(-1i*pi*(ws(k)+sign(ws(k))*(boundaries(l+1) ...
                  -2*boundaries(l)))/(2*wl))/sqrt(wl);
           end
        end
    end
    
    % last one
    mfb{end}=zeros(length(ws),1);

    wl=pi-boundaries(end);

    for k=1:length(mfb{end})
       if boundaries(end)<=abs(ws(k))
             mfb{end}(k)=exp(-1i*pi*(ws(k)+sign(ws(k))*(pi ...
                 -2*boundaries(end)))/(2*wl))/sqrt(wl);
       end
    end
else
    boundaries(boundaries>=pi)=-2*pi+boundaries(boundaries>=pi);
    % note that, contrarily to the papers, the negative boundaries are
    % listed AFTER the positive ones.
    
    mfb=cell(Npic,1);

    % scaling function
    mfb{1}=zeros(length(ws),1);

    for k=1:length(mfb{1})
       if (boundaries(end)<=ws(k)) && (ws(k)<boundaries(1))
          mfb{1}(k)=1;
       end
    end
    
    % each wavelet (High frequencies in a separate way)
    shift=1;
    for l=1:Npic-1
        if boundaries(l)>boundaries(l+1) %high pass filters
            %positive frequencies
            mfb{l+shift}=zeros(length(ws),1);
            wl=pi-boundaries(l);

            for k=1:length(mfb{l+shift})
               if boundaries(l)<=ws(k)
                   mfb{l+shift}(k)=exp(-1i*pi*(ws(k)+sign(ws(k))*(pi ...
                       -2*boundaries(l)))/(2*wl))/sqrt(wl);
               end
            end
            shift=2;
            %negative frequencies
            mfb{l+shift}=zeros(length(ws),1);
            wl=boundaries(l+1)+pi;

            for k=1:length(mfb{l+shift})
               if ws(k)<boundaries(l+1)
                   mfb{l+shift}(k)=exp(-1i*pi*(ws(k)-2*boundaries(l+1) ...
                       -pi)/(2*wl))/sqrt(wl);
               end
            end
        else %individual wavelets
           mfb{l+shift}=zeros(length(ws),1);

            wl=boundaries(l+1)-boundaries(l);

            for k=1:length(mfb{l+shift})
                if (boundaries(l)<=ws(k)) && (ws(k)<boundaries(l+1))
                    if boundaries(l)>0
                        mfb{l+shift}(k)=exp(-1i*pi*(ws(k)+(boundaries(l+1) ...
                            -2*boundaries(l)))/(2*wl))/sqrt(wl);
                    else
                        mfb{l+shift}(k)=exp(-1i*pi*(ws(k)-2*boundaries(l+1) ...
                            +boundaries(l))/(2*wl))/sqrt(wl);
                    end
                end
            end
        end
    end
end