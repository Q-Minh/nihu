function [Tnear, Tfar] = build_block_tree(Ctree)
%BUILD_BLOCK_TREE build block structure from a cluster tree
%   T = BUILD_BLOCK_TREE(C) builds the block structure of the ACA method
%   from the cluster tree C. T is matrix with two columns, each row
%   represents a pair of clusters that needs to be taken into account in
%   the low rank approximation of the matrix.
%
%Example:
%  m = create_sphere_boundary(1, 15);
%  Ctree = build_cluster_tree(m.Nodes(:,2:4));
%  Btree = build_block_tree(Ctree);
%
% See also: build_cluster_tree
%
% Copyright (C) 2014 Peter Fiala

% preallocating the output
CapacityFar = 5000;
Tfar = zeros(CapacityFar,2);
EndFar = 0;

CapacityNear = 5000;
Tnear = zeros(CapacityNear,2);
EndNear= 0;

block_partition(1,1);

Tfar = Tfar(1:EndFar,:);
Tnear = Tnear(1:EndNear,:);

    function block_partition(i,j)
        sigma = Ctree(i);
        tau = Ctree(j);
        if is_admissible(sigma, tau)
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
    end

end

function b = is_admissible(C1, C2, eta)
if nargin < 3
    eta = .5;
end
bb1 = C1.bb;
bb2 = C2.bb;
dist = sqrt(...
    min(min(abs(bsxfun(@minus, bb1(:,1), bb2(:,1)'))))^2 +...
    min(min(abs(bsxfun(@minus, bb1(:,2), bb2(:,2)'))))^2 +...
    min(min(abs(bsxfun(@minus, bb1(:,3), bb2(:,3)'))))^2);
b = max(C1.D, C2.D) < eta * dist;
end
