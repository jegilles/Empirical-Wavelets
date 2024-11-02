function [mapping,risks] = EWT2D_demons_mappings(Lambda,supports,demons)

%==========================================================================
% function [Omega,phi,phi_Omega,interp_function] = EWT2D_demons_mappings(Lambda,supports,demons)
%
% This function returns a partition of the Fourier domain following the
% Voronoi or Watershed method.
%
% Input:
%   -Lambda: mother wavelet's Fourier support
%   -supports: a label image corresponding to the partitioning of the spectrum
%   -demons: structure containing the demons parameters
%       - 'algo': 'Thirion' or 'Vercauteren' for the type of demons
%       - 'level': number of levels in the multiresolution Demons scheme
%       - 'diffus': diffusion regularization parameter in the demons
%       - 'noise': noise parameter in the Vercautern's demons
%       - 'update': 'additive' or 'compositive' for the Vercauteren's demons update
%       - 'field': 'identity' or 'exponential' for the field in the Vercauteren's demons
%
% Output:
%       - Omega: Fourier supports
%       - phi: estimated mappings
%       - phi_Omega: mapped set Lambda
%       - interp_function: 'imwarp or 'iminterpolate'
%
% Author: Charles-Gerard Lucas
% Institution: SDSU - Department of Mathematics
% Year: 2024
% Version: 1.0

if ~isfield(demons,'algo'), demons.algo = 'Vercauteren'; end
if ~isfield(demons,'level'), demons.level = max_level-1:max_level; end
if ~isfield(demons,'diffus'), demons.diffus = .3:.01:.5; end
if ~isfield(demons,'field'), demons.field = 'identity'; end

% Fourier supports
num_supports = max(supports,[],'all');
mapping.Omega = cell(1,num_supports);
for k=1:num_supports, mapping.Omega{k} = double(bwlabel((supports==k),8)==1); end

% Demons algorithm
switch demons.algo
    case 'Thirion'
        demons.interp_type = 'linear';
        demons.interp_function = 'imwarp';
    case 'Vercauteren'
        if ~isfield(demons,'noise'), demons.noise = 4; end
        if ~isfield(demons,'update'), demons.update = 'additive'; end
        if strcmp(demons.field,'exp') 
            demons.interp_type='spline'; 
        else 
            demons.interp_type='linear';
        end
        demons.interp_function = 'iminterpolate';
end
mapping.interp_function = demons.interp_function;
mapping.interp_type = demons.interp_type;

% Grid search on demons
mapping.phi = cell(1,num_supports); mapping.phi_Omega = cell(1,num_supports); risks = cell(1,num_supports);
for k=1:num_supports, [mapping.phi{k},mapping.phi_Omega{k},risks{k}] = EWT2D_demons_grid_search(Lambda,mapping.Omega{k},demons); end