function Show_EWT_Boundaries(f,boundaries,R,SamplingRate,logspec,InitBounds,presig)

%===============================================================================
% function Show_EWT_Boundaries(f,boundaries,R,SamplingRate,logspec,InitBounds,presig)
% 
% This function plot the boundaries by superposition on the graph 
% of the magnitude of the Fourier spectrum on the frequency interval 
% ([0,pi] for a real signal and [0,2pi] for a complex signal. 
% If the sampling rate is provided, then the horizontal axis
% will correspond to the real frequencies. If provided, the plot will
% superimposed a set of initial boundaries (in black). The input presig 
% will be superimposed on the original plot (useful to visualize regularized 
% version of f)
%
% Input:
%   - f: signal
%   - boundaries: list of all boundaries
%   - R: ratio to plot only a portion of the spectrum (we plot the
%        interval [0,pi/R]. R must be >=1. This applies for the real signal
%        case only.
%   - SamplingRate: sampling rate used during the signal acquisition (must
%                   be set to -1 if it is unknown)
%		- InitBounds: initial bounds when you use an adaptive detection method
%		- presig: preprocessed version of the spectrum on which the detection is made
%   - logspec: flag to indicate if the logarithm of the spectrum must be
%              displayed (set the flag to 1, any other value will plot the regular
%              spectrum)
%
% Optional inputs:
%   - InitBounds: optional initial boundaries
%   - presig: preprocessed signal
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 2.0
%===============================================================================

figure;
freq=2*pi*[0:length(f)-1]/length(f);

if SamplingRate~=-1
    freq=freq*SamplingRate/(2*pi);
    boundaries=boundaries*SamplingRate/(2*pi);
end

if isreal(f)
    if R<1
        R=1;
    end
    %plot only half the spectrum in the case
    R=round(length(f)/(2*R));
else
    R=length(f);
end

if logspec==1
        f=log(abs(fft(f)));
    else
        f=abs(fft(f));
end

plot(freq(1:R),f(1:R));
hold on
if nargin>6
   plot(freq(1:R),presig(1:R),'color',[0,0.5,0],'LineWidth',2); 
end
NbBound=length(boundaries);
 
for i=1:NbBound
     if boundaries(i)>freq(R)
         break
     end
     line([boundaries(i) boundaries(i)],[0 max(f)],'LineWidth',2,'LineStyle','--','Color',[1 0 0]);
end
 
if nargin>5
    NbBound=length(InitBounds);
    InitBounds=InitBounds*2*pi/length(f);
    
    for i=1:NbBound
        if InitBounds(i)>freq(R)
            break
        end
        line([InitBounds(i) InitBounds(i)],[0 max(f)],'LineWidth',1,'LineStyle','--','Color',[0 0 0]);
    end
end