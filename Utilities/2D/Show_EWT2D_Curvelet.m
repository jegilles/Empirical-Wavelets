function Show_EWT2D_Curvelet(ewtc,option)

%==========================================================================
% function Show_EWT2D_Curvelet(ewtc,option)
% 
% This function displays the curvelet coefficient obtained by the Empirical 
% Curvelet Transform (one figure per scale)
%
% Inputs:
%   -ewtc: output of the empirical curvelet transform
%   -option: type of curvelet transform (see EWT2D_Curvelet documentation)
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics & Statistics
% Version: 1.0 (2013)
% Version: 2.0 (2019)
%==========================================================================

Ns=length(ewtc);

% Show the low pass filter
figure;
imshow(ewtc{1},[]);title('s=1');

if option<3
    for s=2:Ns
       Nt=length(ewtc{s});
       figure;
       Nr=floor(sqrt(Nt));
       Nc=Nt-Nr^2+Nr;
       for t=1:Nt
          subplot(Nr,Nc,t);
          imshow(ewtc{s}{t},[]); title(['s = ' num2str(s) ', t = ' num2str(t)]);
       end
    end
else
    for s=2:Ns
       Nt=length(ewtc{s});
       figure;
       Nr=floor(sqrt(Nt));
       Nc=Nt-Nr^2+Nr;
       for t=1:Nt
          subplot(Nr,Nc,t);
          imshow(ewtc{s}{t},[]);title(['s = ' num2str(t) ', t = ' num2str(s)]);
       end
    end
end
