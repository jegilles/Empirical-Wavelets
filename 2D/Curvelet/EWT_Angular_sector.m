function Angular = EWT_Angular_sector(theta,radius,theta0,theta1,r0,r1,gammaw,Dtheta)

%==================================================================================
% function Angular = EWT_Angular_sector(theta,radius,theta0,theta1,r0,r1,gammaw,Dtheta)
% 
% This function creates the curvelet filter in the Fourier domain which
% have a support defined in polar coordinates (r,angles) in
% [(1-gammaw)r0,(1+gammaw)r1]x[theta0-Dtheta,theta1+Dtheta] and 
% [(1-gammaw)r0,(1+gammaw)r1]x[theta0-Dtheta+pi,theta1+Dtheta+pi]
%
% Inputs:
%   -theta: angle grid
%   -radius: radius grid
%   -theta0,theta1: the consecutive angles defining the support's limits
%   (theta0 < theta1)
%   -r0,r1: the consecutive scales defining the support's limits
%   (r0<r1<=pi)
%   -gammaw: coefficient defining the spread of scale transition zones
%   -Dtheta: coefficient defining the spread of angle transition zones
%
% Output:
%   -Angular: curvelet filter in Fourier domain
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics and Statistics
% Version: 1.0 (2013)
% Version: 2.0 (2015)
% Version: 3.0 (2019)
%==================================================================================

H=size(theta,2);
W=size(theta,1);

wan=1/(2*gammaw*r0);
wam=1/(2*gammaw*r1);
wpbn=(1+gammaw)*r0;
wmbn=(1-gammaw)*r0;
wpbm=(1+gammaw)*r1;
wmbm=(1-gammaw)*r1;

an=1/(2*Dtheta);
pbn=theta0+Dtheta;
mbn=theta0-Dtheta;
pbm=theta1+Dtheta;
mbm=theta1-Dtheta;


Angular=zeros(size(theta));

if r1<pi
    for i=1:W
       for j=1:H
            if ((theta(j,i)>pbn) && (theta(j,i)<mbm))
                if ((radius(j,i)>wpbn) && (radius(j,i)<wmbm)) %inside
                    Angular(j,i)=1;
                elseif ((radius(j,i)>=wmbm) && (radius(j,i)<=wpbm)) %top of inside - radial only
                    Angular(j,i)=cos(pi*EWT_beta(wam*(radius(j,i)-wmbm))/2);
                elseif ((radius(j,i)>=wmbn) && (radius(j,i)<=wpbn)) %bottom of inside - radial only
                    Angular(j,i)=sin(pi*EWT_beta(wan*(radius(j,i)-wmbn))/2);
                end
            elseif ((theta(j,i)>=mbn) && (theta(j,i)<=pbn)) 
                if ((radius(j,i)>wpbn) && (radius(j,i)<wmbm)) %left of inside - angular only
                    Angular(j,i)=sin(pi*EWT_beta(an*(theta(j,i)-mbn))/2);
                elseif ((radius(j,i)>wmbm) && (radius(j,i)<=wpbm)) %top-left of inside - radial/angular mix
                    Angular(j,i)=sin(pi*EWT_beta(an*(theta(j,i)-mbn))/2)*cos(pi*EWT_beta(wam*(radius(j,i)-wmbm))/2);
                elseif ((radius(j,i)>=wmbn) && (radius(j,i)<=wpbn)) %bottom-left of inside - radial/angular mix
                    Angular(j,i)=sin(pi*EWT_beta(an*(theta(j,i)-mbn))/2)*sin(pi*EWT_beta(wan*(radius(j,i)-wmbn))/2);
                end
            elseif ((theta(j,i)>=mbm) && (theta(j,i)<=pbm))
                if ((radius(j,i)>wpbn) && (radius(j,i)<wmbm)) %right of inside - angular only
                    Angular(j,i)=cos(pi*EWT_beta(an*(theta(j,i)-mbm))/2);
                elseif ((radius(j,i)>wmbm) && (radius(j,i)<=wpbm)) %top-right of inside - radial/angular mix
                    Angular(j,i)=cos(pi*EWT_beta(an*(theta(j,i)-mbm))/2)*cos(pi*EWT_beta(wam*(radius(j,i)-wmbm))/2);
                elseif ((radius(j,i)>=wmbn) && (radius(j,i)<=wpbn)) %bottom-right of inside - radial/angular mix
                    Angular(j,i)=cos(pi*EWT_beta(an*(theta(j,i)-mbm))/2)*sin(pi*EWT_beta(wan*(radius(j,i)-wmbn))/2);
                end
            end
       end
   end
else
   for i=1:W
       for j=1:H
            if ((theta(j,i)>pbn) && (theta(j,i)<mbm))
                if (radius(j,i)>wpbn) %inside
                    Angular(j,i)=1;
                elseif ((radius(j,i)>=wmbn) && (radius(j,i)<=wpbn)) %bottom of inside - radial only
                    Angular(j,i)=sin(pi*EWT_beta(wan*(radius(j,i)-wmbn))/2);
                end
             elseif ((theta(j,i)>=mbn) && (theta(j,i)<=pbn))
                if (radius(j,i)>wpbn) %left of inside - angular only
                    Angular(j,i)=sin(pi*EWT_beta(an*(theta(j,i)-mbn))/2);
                elseif ((radius(j,i)>=wmbn) && (radius(j,i)<=wpbn)) %bottom-left of inside - radial/angular mix
                    Angular(j,i)=sin(pi*EWT_beta(an*(theta(j,i)-mbn))/2)*sin(pi*EWT_beta(wan*(radius(j,i)-wmbn))/2);
                end
            elseif ((theta(j,i)>=mbm) && (theta(j,i)<=pbm)) 
                if (radius(j,i)>wpbn) %right of inside - angular only
                    Angular(j,i)=cos(pi*EWT_beta(an*(theta(j,i)-mbm))/2);
                elseif ((radius(j,i)>=wmbn) && (radius(j,i)<=wpbn)) %bottom-right of inside - radial/angular mix
                    Angular(j,i)=cos(pi*EWT_beta(an*(theta(j,i)-mbm))/2)*sin(pi*EWT_beta(wan*(radius(j,i)-wmbn))/2);
                end
            end
       end
   end
end
% we symmetrize the filter
Angular = Angular + EWT2D_Symmterize_Curvelet_Filter(Angular);