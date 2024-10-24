function ouv=EWT_FunctionOpening(f,sizeel)

%=========================================================
% function ouv=EWT_FunctionOpening(f,sizeel)
%
% This function perform the mathematical morphology 
% opening (= dilation of erosion of f) operator for 
% functions according to structural element of size sizeel.
%
% Inputs:
%   -f: input function
%   -sizeel: size of the structural element
%
% Output:
%   -ouv: the opened function
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
%=========================================================
ouv=EWT_FunctionDilation(EWT_FunctionErosion(f,sizeel),sizeel);