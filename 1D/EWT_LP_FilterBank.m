function mfb=EWT_LP_FilterBank(boundaries,N,cpx)

% =========================================================================
% function mfb=EWT_LP_FilterBank(boundaries,N,cpx)
%
% This function generate the Littlewood-Paley filter bank (scaling function
% + wavelets) corresponding to the provided set of frequency segments
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
% Year: 2019
% Version: 2.0
% =========================================================================

Npic=length(boundaries);

%create the frequency axis [0,pi)u[-pi,0)
Mi=floor(N/2)+1;
w=(0:2*pi/N:2*pi-2*pi/N)';
w(Mi+1:end)=-2*pi+w(Mi+1:end); 

if cpx==0
    % We compute gamma accordingly to the theory
    r=[1;([boundaries(2:end) ; pi]- boundaries)./([boundaries(2:end) ; pi]+boundaries)];
    gamma=min(r);
    gamma=(1-1/N)*gamma; %this ensure that gamma is chosen as strictly less than the min
    
    %make the second half of the axis be negative frequencies
    mfb=cell(Npic+1,1);
    % We start by generating the scaling function
    mfb{1}=zeros(N,1);
    
    am=1/(2*gamma*abs(boundaries(1)));
    pbm=(1+gamma)*boundaries(1);
    mbm=(1-gamma)*boundaries(1);
    for k=1:N
       if abs(w(k))<=mbm
           mfb{1}(k)=1;
       elseif ((abs(w(k))>=mbm) && (abs(w(k))<=pbm))
           mfb{1}(k)=cos(pi*EWT_beta(am*(abs(w(k))-mbm))/2);
       end
    end

    % We generate the wavelets except the last one
    for l=1:Npic-1
        mfb{l+1}=zeros(N,1);
        an=1/(2*gamma*abs(boundaries(l)));
        pbn=(1+gamma)*boundaries(l);
        mbn=(1-gamma)*boundaries(l);

        am=1/(2*gamma*abs(boundaries(l+1)));
        pbm=(1+gamma)*boundaries(l+1);
        mbm=(1-gamma)*boundaries(l+1);

        for k=1:N
           if ((abs(w(k))<=mbm) && (abs(w(k))>=pbn))
               mfb{l+1}(k)=1;
           elseif ((abs(w(k))>=mbm) && (abs(w(k))<=pbm))
               mfb{l+1}(k)=cos(pi*EWT_beta(am*(abs(w(k))-mbm))/2);
           elseif ((abs(w(k))>=mbn) && (abs(w(k))<=pbn))
               mfb{l+1}(k)=sin(pi*EWT_beta(an*(abs(w(k))-mbn))/2);
           end
        end
    end
    
    % last one
    mfb{Npic+1}=zeros(N,1);
    am=1/(2*gamma*abs(boundaries(end)));
    mbm=(1-gamma)*boundaries(end);
    pbm=(1+gamma)*boundaries(end);
    for k=1:N
       if abs(w(k))>=pbm
           mfb{Npic+1}(k)=1;
       elseif ((abs(w(k))>=mbm) && (abs(w(k))<=pbm))
           mfb{Npic+1}(k)=sin(pi*EWT_beta(am*(abs(w(k))-mbm))/2);
       end
    end
else
    boundaries(boundaries>=pi)=-2*pi+boundaries(boundaries>=pi);
    sb=[boundaries(boundaries<0) ; boundaries(boundaries>=0)];
    r=([sb ; pi]-[-pi ; sb])./abs([sb ; pi]+[-pi ; sb]);
    gamma=min(r);
    gamma=(1-1/N)*gamma; %this ensure that gamma is chosen as strictly less than the min

    
    mfb=cell(Npic+1,1);
    % We start by generating the scaling function 
    %(here we assume the 0 is not a possible boundary - COULD CHANGE IN THE FUTURE)
    mfb{1}=zeros(N,1);
    %compute the bounds of the transition area
    an=1/(2*gamma*abs(boundaries(end)));
    pbn=(1+sign(boundaries(end))*gamma)*boundaries(end);
    mbn=(1-sign(boundaries(end))*gamma)*boundaries(end);
    
    am=1/(2*gamma*abs(boundaries(1)));
    pbm=(1+sign(boundaries(1))*gamma)*boundaries(1);
    mbm=(1-sign(boundaries(1))*gamma)*boundaries(1);

    for k=1:N
       if ((w(k)<=mbm) && (w(k)>=pbn))
           mfb{1}(k)=1;
       elseif ((w(k)>=mbm) && (w(k)<=pbm))
           mfb{1}(k)=cos(pi*EWT_beta(am*(w(k)-mbm))/2);
       elseif ((w(k)>=mbn) && (w(k)<=pbn))
           mfb{1}(k)=sin(pi*EWT_beta(an*(w(k)-mbn))/2);
       end
    end
    
    % We generate the wavelets
    shift=1;
    for l=1:Npic-1
        if boundaries(l)>boundaries(l+1) %high pass filters
            mfb{l+shift}=zeros(N,1);
            %compute the bounds of the transition area
            an=1/(2*gamma*abs(boundaries(l)));
            pbn=(1+sign(boundaries(l))*gamma)*boundaries(l);
            mbn=(1-sign(boundaries(l))*gamma)*boundaries(l);
            for k=1:N
               if w(k)>=pbn
                   mfb{l+shift}(k)=1;
               elseif ((w(k)>=mbn) && (w(k)<=pbn))
                   mfb{l+shift}(k)=sin(pi*EWT_beta(an*(w(k)-mbn))/2);
               end
            end
            
            shift=2;
            mfb{l+shift}=zeros(N,1);
            %compute the bounds of the transition area
            am=1/(2*gamma*abs(boundaries(l+1)));
            pbm=(1+sign(boundaries(l+1))*gamma)*boundaries(l+1);
            mbm=(1-sign(boundaries(l+1))*gamma)*boundaries(l+1);
            for k=1:N
               if w(k)<=mbm
                   mfb{l+shift}(k)=1;
               elseif ((w(k)>=mbm) && (w(k)<=pbm))
                   mfb{l+shift}(k)=cos(pi*EWT_beta(am*(w(k)-mbm))/2);
               end
            end            
        else 
            mfb{l+shift}=zeros(N,1);
            %compute the bounds of the transition area
            an=1/(2*gamma*abs(boundaries(l)));
            pbn=(1+sign(boundaries(l))*gamma)*boundaries(l);
            mbn=(1-sign(boundaries(l))*gamma)*boundaries(l);


            am=1/(2*gamma*abs(boundaries(l+1)));
            pbm=(1+sign(boundaries(l+1))*gamma)*boundaries(l+1);
            mbm=(1-sign(boundaries(l+1))*gamma)*boundaries(l+1);

            for k=1:N
               if ((w(k)<=mbm) && (w(k)>=pbn))
                   mfb{l+shift}(k)=1;
               elseif ((w(k)>=mbm) && (w(k)<=pbm))
                   mfb{l+shift}(k)=cos(pi*EWT_beta(am*(w(k)-mbm))/2);
               elseif ((w(k)>=mbn) && (w(k)<=pbn))
                   mfb{l+shift}(k)=sin(pi*EWT_beta(an*(w(k)-mbn))/2);
               end
            end
        end
    end
end