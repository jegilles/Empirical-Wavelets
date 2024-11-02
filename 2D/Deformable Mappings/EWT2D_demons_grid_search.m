function [phi,phi_Omega,risks] = EWT2D_demons_grid_search(Lambda,Omega,demons)

%==========================================================================
% This function performs a grid search on the demons algorithm.
%   
% Input:
%   -Lambda: Moving support
%   -Omega: Fixed support
%   -demons: structure containing
%        - 'algo': 'Thirion' or 'Vercauteren' for the type of Demons
%        - 'level': number of levels in the multiresolution Demons scheme
%        - 'diffus': diffusion regularization parameter in the Demons
%        - 'spatial_uncertainty': spatial uncertainty parameter in the Vercautern's Demons
%        - 'update': 'additive' or 'compositive' for the Vercauteren's Demons update
%        - 'field': 'identity' or 'exponential' for the field in the
%        Vercauteren's Demons
%       - 'interp_type': 'identity' or 'exponential' for the field in the
%
% Output:
%   -mfb: constructed filterbank
%
% Author: Charles-Gerard Lucas
% Institution: SDSU - Department of Mathematics
% Year: 2024
% Version: 1.0
%==========================================================================

switch demons.algo
    case 'Thirion'
        s1 = length(demons.diffus); s2 = length(demons.level);
        phi_tmp = cell(s1,s2); phi_Omega_tmp = cell(s1,s2);
        risks = zeros(s1,s2);
        for l1=1:length(demons.diffus)
            for l2=1:length(demons.level)
                niter = 2^3*2.^(demons.level(l2):-1:1);
                [phi_tmp{l1,l2},phi_Omega_tmp{l1,l2}] = imregdemons(Lambda,Omega,niter,'AccumulatedFieldSmoothing',...
                    demons.diffus(l1),'PyramidLevels',demons.level(l2),'DisplayWaitbar',false);
                risks(l1,l2) = norm(phi_Omega_tmp{l1,l2}-Omega,'fro');
            end
        end
        [la,lp] = find(risks==min(risks,[],'all'),1);
        phi = phi_tmp{la,lp};
        phi_Omega = phi_Omega_tmp{la,lp};
    case 'Vercauteren'
        s1 = length(demons.diffus); s2 = length(demons.spatial_uncertainty); s3 = length(demons.level);
        phi_tmp = cell([s1,s2,s3]); phi_Omega_tmp = cell([s1,s2,s3]);
        risks = zeros([s1,s2,s3]);
        for l1=1:s1
            for l2=1:s2
                for l3=1:s3
                    params = struct();
                    params.sigma_i = 1;
                    params.sigma_fluid = 1; 
                    params.sigma_diffusion = demons.diffus(l1);
                    params.sigma_x = demons.spatial_uncertainty(l2);
                    params.nlevel = demons.level(l3);
                    params.do_display = 0;
                    params.interp_type = demons.interp_type;
                    params.update = demons.update;
                    params.field = demons.field;
                    params.niter = 500;
                    params.stop_criterium = 1e-3;
                    [phi_tmp{l1,l2,l3},phi_Omega_tmp{l1,l2,l3}] = Vercauteren_demons(Lambda,Omega,params);
                    risks(l1,l2,l3) = norm(phi_Omega_tmp{l1,l2,l3}./256-Omega,'fro')^2;
                end
            end
        end
        [~,loc] = min(risks(:)); [ld,ln,ll] = ind2sub(size(risks),loc);
        phi = phi_tmp{ld,ln,ll}; phi_Omega = phi_Omega_tmp{ld,ln,ll};
end

end