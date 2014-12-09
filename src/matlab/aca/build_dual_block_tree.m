function [Tnear, Tfar] = build_dual_block_tree(RowTree, ColTree, is_admissible)
%BUILD_DUAL_BLOCK_TREE build block structure from two cluster trees
%   [Tnear, Tfar] = BUILD_DUAL_BLOCK_TREE(RowTree, ColTree, Admit) builds the
%   block structure from the cluster trees RwoTree and ColTree.
%   Admit is a function deciding if two clusters are admissible to form a
%   block or not.
%   Tnear and Tfar are matrices with two columns, each row
%   representing a pair of clusters.
%
%Example:
%  r = create_sphere_boundary(1, 15);
%  c = translate_mesh(r, [2 2 2]);
%  RowTree = build_cluster_tree(r.Nodes(:,2:4));
%  ColTree = build_cluster_tree(c.Nodes(:,2:4));
%  [Bnear, Bfar] = build_block_tree(RowTree, ColTree);
%
% See also: build_cluster_tree build_block_tree

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
        sigma = RowTree(i);
        tau = ColTree(j);
        if is_admissible(sigma, tau)
            EndFar = EndFar+1;
            if EndFar > CapacityFar
                CapacityFar = 2*CapacityFar;
                Tfar(CapacityFar,2) = 0;
            end
            Tfar(EndFar,:) = [i, j];
        elseif isempty(sigma.children) || isempty(tau.children)
            EndNear = EndNear+1;
            if EndNear > CapacityNear
                CapacityNear = 2*CapacityNear;
                Tnear(CapacityNear,2) = 0;
            end
            Tnear(EndNear,:) = [i, j];
        else
            [s, t] = divide_both(i, j);
            for i = 1 : numel(s)
                for j = 1 : numel(t)
                    block_partition(s(i), t(j));
                end
            end
        end % if
    end % of function block_partition

    function [s, t] = divide_larger(i, j)
        sigma = RowTree(i);
        tau = ColTree(j);
        % try to divide the larger cluster if possible
        Ds = norm(diff(sigma.bb));	% diameter of sigma
        Dt = norm(diff(tau.bb));	% diameter of tau
        if (Ds > Dt && ~isempty(sigma.children)) || isempty(tau.children) % divide sigma
            s = sigma.children;
            t = j;
        else % divide tau
            s = i;
            t = tau.children;
        end
    end

    function [s, t] = divide_both(i, j)
        s = RowTree(i).children;
        t = ColTree(j).children;
    end
end % of function build_dual_block_tree
