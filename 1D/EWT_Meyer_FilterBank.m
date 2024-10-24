function [mfb,ws]=EWT_Meyer_FilterBank(boundaries,N,cpx)

% =========================================================================
% function mfb=EWT_Meyer_FilterBank(boundaries,N,cpx)
%
% This function generate the Meyer filter bank (scaling function
% + wavelets) corresponding to the provided set of frequency segments.
% In this version, the highpass filters are made of half a Meyer filter and
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


if cpx == 0 %%REAL CASE
    mfb=cell(Npic+1,1);

    % scaling function
    mfb{1}=zeros(length(ws),1);
    wm=(boundaries(1)+boundaries(2))/2;
    for k=1:length(mfb{1})
        if abs(ws(k))<=wm
            mfb{1}(k)=cos(pi*EWT_beta(abs(ws(k))/wm)/2)/sqrt(wm);
        end
    end
    
    % each wavelet except the last one
    for l=1:Npic-1
       mfb{l+1}=zeros(length(ws),1);

        wn=(boundaries(l+1)+boundaries(l))/2;
        if l==1
            wm=boundaries(l)/2;
        else
            wm=(boundaries(l-1)+boundaries(l))/2;
        end
        if l==(Npic-1)
            wp=(pi+boundaries(l+1))/2;
        else
            wp=(boundaries(l+1)+boundaries(l+2))/2;
        end
        for k=1:length(mfb{l+1})
            if (abs(ws(k))>=wm) && (abs(ws(k))<=wn)
                mfb{l+1}(k)=sqrt(2/(wp-wm))*exp(1i*pi*sign(ws(k))*(abs(ws(k))+(wp-4*wm)/3)/(wp-wm))*sin(pi*EWT_beta((abs(ws(k))-wm)/(wn-wm))/2);
            elseif (abs(ws(k))>=wn) && (abs(ws(k))<=wp)
                mfb{l+1}(k)=sqrt(2/(wp-wm))*exp(1i*pi*sign(ws(k))*(abs(ws(k))+(wp-4*wm)/3)/(wp-wm))*cos(pi*EWT_beta((abs(ws(k))-wn)/(wp-wn))/2);                
            end
        end
    end
    
    % last one
    mfb{end}=zeros(length(ws),1);

    wn=(pi+boundaries(end))/2;
    wm=(boundaries(end)+boundaries(end-1))/2;

    for k=1:length(mfb{end})
        if (wm<=abs(ws(k))) && (abs(ws(k))<=wn)
            mfb{end}(k)=sqrt(2/(wn-wm))*exp(1i*pi*sign(ws(k))*(abs(ws(k))+(pi-4*wm)/3)/(pi-wm))*sin(pi*EWT_beta((abs(ws(k))-wm)/(wn-wm))/2);
        elseif abs(ws(k))>=wn
            mfb{end}(k)=sqrt(2/(wn-wm))*exp(1i*pi*sign(ws(k))*(abs(ws(k))+(pi-4*wm)/3)/(pi-wm));
        end
    end
else %%COMPLEX CASE
    boundaries(boundaries>=pi)=-2*pi+boundaries(boundaries>=pi);
    % note that, contrarily to the papers, the negative boundaries are
    % listed AFTER the positive ones.
    
    mfb=cell(Npic,1);

    % scaling function
    mfb{1}=zeros(length(ws),1);
    
    wn=(boundaries(end)+boundaries(1))/2;
    wm=(boundaries(end-1)+boundaries(end))/2;
    wp=(boundaries(1)+boundaries(2))/2;
    for k=1:length(mfb{1})
        if (ws(k)>=wm) && (ws(k)<=wn)
            mfb{1}(k)=sqrt(2/(wp-wm))*sin(pi*EWT_beta((ws(k)-wm)/(wn-wm))/2);
        elseif (ws(k)>=wn) && (ws(k)<=wp)
            mfb{1}(k)=sqrt(2/(wp-wm))*cos(pi*EWT_beta((ws(k)-wn)/(wp-wn))/2);
        end
    end
    
    % each wavelet (High frequencies in a separate way)
    shift=1;
    for l=1:Npic-1
        if boundaries(l)>boundaries(l+1) %high pass filters
            %positive frequencies
            mfb{l+shift}=zeros(length(ws),1);

            wn=(pi+boundaries(l))/2;
            if l==1
                wm=(boundaries(end)+boundaries(l))/2;
            else
                wm=(boundaries(l)+boundaries(l-1))/2;
            end

            for k=1:length(mfb{l+shift})
                if (wm<=ws(k)) && (ws(k)<=wn)
                    if ws(k)>=0
                        mfb{l+shift}(k)=sqrt(2/(wn-wm))*exp(1i*pi*(ws(k)+(pi-4*wm)/3)/(pi-wm))*sin(pi*EWT_beta((ws(k)-wm)/(wn-wm))/2);
                    else
                        mfb{l+shift}(k)=sqrt(2/(wn-wm))*exp(1i*pi*(ws(k)+(wm-4*pi)/3)/(pi-wm))*sin(pi*EWT_beta((ws(k)-wm)/(wn-wm))/2);
                    end
                elseif ws(k)>=wn
                    if ws(k)>=0
                        mfb{l+shift}(k)=sqrt(2/(wn-wm))*exp(1i*pi*(ws(k)+(pi-4*wm)/3)/(pi-wm));
                    else
                        mfb{l+shift}(k)=sqrt(2/(wn-wm))*exp(1i*pi*(ws(k)+(wm-4*pi)/3)/(pi-wm));
                    end
                end
            end
            shift=2;
            %negative frequencies
            mfb{l+shift}=zeros(length(ws),1);
            
            wn=(-pi+boundaries(l+1))/2;
            if (l+1)==Npic
                wp=(boundaries(1)+boundaries(l+1))/2;
            else
                wp=(boundaries(l+1)+boundaries(l+2))/2;
            end

            for k=1:length(mfb{l+shift})
                if (wn<=ws(k)) && (ws(k)<=wp)
                    if ws(k)>=0
                        mfb{l+shift}(k)=sqrt(2/(wp-wn))*exp(1i*pi*(ws(k)+(wp+4*pi)/3)/(wp+pi))*cos(pi*EWT_beta((ws(k)-wn)/(wp-wn))/2);
                    else
                        mfb{l+shift}(k)=sqrt(2/(wp-wn))*exp(1i*pi*(ws(k)-(4*wp+pi)/3)/(wp+pi))*cos(pi*EWT_beta((ws(k)-wn)/(wp-wn))/2);
                    end
                elseif ws(k)<=wn
                    if ws(k)>=0
                        mfb{l+shift}(k)=sqrt(2/(wp-wn))*exp(1i*pi*(ws(k)+(wp+4*pi)/3)/(wp+pi));
                    else
                        mfb{l+shift}(k)=sqrt(2/(wp-wn))*exp(1i*pi*(ws(k)-(4*wp+pi)/3)/(wp+pi));
                    end
                end
            end
        else %individual wavelets
            mfb{l+shift}=zeros(length(ws),1);

            wn=(boundaries(l+1)+boundaries(l))/2;
            
            if (l-1)<=0
                wm=(boundaries(l)+boundaries(end))/2;
            else
                if (boundaries(l-1)>0) && (boundaries(l)<0)
                    wm=(boundaries(l)-pi)/2;
                else
                    wm=(boundaries(l)+boundaries(l-1))/2;
                end
            end
            
            if (l+2)>=Npic
                wp=(boundaries(l+1)+boundaries(1))/2;
            else
                if (boundaries(l+1)>0) && (boundaries(l+2)<0)
                    wp=(boundaries(l+1)+pi)/2;
                else
                    wp=(boundaries(l+1)+boundaries(l+2))/2;
                end
            end
            
            for k=1:length(mfb{l+1})
                if (ws(k)>=wm) && (ws(k)<=wn)
                    if sign(ws(k))>=0
                        mfb{l+shift}(k)=sqrt(2/(wp-wm))*exp(1i*pi*(ws(k)+(wp-4*wm)/3)/(wp-wm))*sin(pi*EWT_beta((ws(k)-wm)/(wn-wm))/2);
                    else
                        mfb{l+shift}(k)=sqrt(2/(wp-wm))*exp(1i*pi*(ws(k)-(4*wp-wm)/3)/(wp-wm))*sin(pi*EWT_beta((ws(k)-wm)/(wn-wm))/2);
                    end
                elseif (ws(k)>=wn) && (ws(k)<=wp)
                    if sign(ws(k))>=0
                        mfb{l+shift}(k)=sqrt(2/(wp-wm))*exp(1i*pi*(ws(k)+(wp-4*wm)/3)/(wp-wm))*cos(pi*EWT_beta((ws(k)-wn)/(wp-wn))/2);
                    else
                        mfb{l+shift}(k)=sqrt(2/(wp-wm))*exp(1i*pi*(ws(k)-(4*wp-wm)/3)/(wp-wm))*cos(pi*EWT_beta((ws(k)-wn)/(wp-wn))/2);
                    end
                end
            end
        end
    end
end