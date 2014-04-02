function T = build_block_tree(Ctree)
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
Capacity = 5000;
T = zeros(Capacity,2);
End = 0;

% preallocating the temporary pairs
T1 = zeros(Capacity,2);

% initial guess pair
T0 = [1 1];

while ~isempty(T0)
    n = 0;  % restart building the temporary pairs
    for i = 1 : size(T0,1)
        if is_admissible(Ctree(T0(i,:)))
            % build into final pairs
            End = End+1;
            if End > Capacity
                Capacity = 2*Capacity;
                T(Capacity,2) = 0;
            end
            T(End,:) = T0(i,:);
        else
            % try with children
            Ch = get_children(Ctree(T0(i,:)));
            l = size(Ch,1);
            T1(n+(1:l),:) = Ch;
            n = n + l;
        end
    end
    % replace guess with children
    T0 = T1(1:n,:);
end

T = T(1:End,:);
end

function TT = get_children(C)
[L, R] = meshgrid(C(1).children, C(2).children);
TT = [L(:), R(:)];
end

function b = is_admissible(C)
b = C(2).D*sqrt(3) < norm(C(1).c - C(2).c) ||...
    length(C(1).ind) == 1 ||...
    length(C(2).ind) == 1;
end
