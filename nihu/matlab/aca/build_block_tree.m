function [Tnear, Tfar] = build_block_tree(Ctree, Admit)
%BUILD_BLOCK_TREE build block structure from a cluster tree
%   [Tnear, Tfar] = BUILD_BLOCK_TREE(C, ADMIT) builds the block structure
%   from the cluster tree C and the function ADMIT.
%   The function ADMIT decides if two clusters are admissible to form a
%   block.
%   Tnear and Tfar are matrices with two columns, each row
%   representing a pair of clusters.
%
%Example:
%  m = create_sphere_boundary(1, 15);
%  Ctree = build_cluster_tree(m.Nodes(:,2:4));
%  [Bnear, Bfar] = build_block_tree(Ctree);
%
% See also: build_cluster_tree

% Copyright (C) 2014-2014 Peter Fiala

% preallocating the output
CapacityFar = 5000;
Tfar = zeros(CapacityFar,2);
EndFar = 0;

CapacityNear = 5000;
Tnear = zeros(CapacityNear,2);
EndNear= 0;

% recursively partition the tree starting with the (1,1) block
block_partition(1,1);

% truncate the block tree
Tfar = Tfar(1:EndFar,:);
Tnear = Tnear(1:EndNear,:);

    function block_partition(i,j)
        sigma = Ctree(i);
        tau = Ctree(j);
        if Admit(sigma, tau)
            EndFar = EndFar+1;
            if EndFar > CapacityFar
                CapacityFar = 2*CapacityFar;
                Tfar(CapacityFar,2) = 0;
            end
            Tfar(EndFar,:) = [i, j];
        elseif ~isempty(sigma.children) && ~isempty(tau.children)
            [s, t] = meshgrid(sigma.children, tau.children);
            s = s(:); t = t(:);
            for k = 1 : length(s)
                block_partition(s(k),t(k));
            end
        else
            EndNear = EndNear+1;
            if EndNear > CapacityNear
                CapacityNear = 2*CapacityNear;
                Tnear(CapacityNear,2) = 0;
            end
            Tnear(EndNear,:) = [i, j];
        end
    end % of internal function block_partition
end % of function build_block_tree
