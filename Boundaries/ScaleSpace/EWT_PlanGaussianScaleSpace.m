function plane=EWT_PlanGaussianScaleSpace(f)

%==========================================================================
%
% function plane=EWT_PlanGaussianScaleSpace(f)
%
% This function builds a scale-space representation of f via successive
% convolutions with a discrete Gaussian Kernel.
%
% Input:
%   -f: the input function
%
% Output:
%   -plane: the scale-space representation of f (the horizontal axis 
%           corresponds to the scales).
%
% Authors: Jerome Gilles & George Istambouli
% Institution: SDSU - Department of Mathematics
% Year: 2016
% Version: 2.0
% =========================================================================
n=3;%4
t=0.16;%initial scale
Niter=ceil(length(f)/n);
ker=besseli(-n:n,t,1);%initialize the discrete Gauusian kernel
plane=sparse(zeros(length(f),Niter+1));%initialize the scale space plane

%find the initial local minima
bounds=EWT_LocalMaxMin2(f,-1);
plane(bounds)=1;
N=zeros(1,Niter+1);
N(1)=length(bounds);

%go through the scale-space
f=f';
for i=1:Niter    
    f=conv(f,ker,'same');%filter   
    bounds=EWT_LocalMaxMin2(f,-1); %local minima
    plane(bounds,i+1)=1;
    N(i+1)=length(bounds);
    if N(i+1)==0
        break
    end
end