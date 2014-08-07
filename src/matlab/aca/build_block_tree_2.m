function [Tnear, Tfar] = build_block_tree_2(SourceTree, ReceiverTree)

% preallocating the output
CapacityFar = 5000;
Tfar = zeros(CapacityFar,2);
EndFar = 0;

CapacityNear = 5000;
Tnear = zeros(CapacityNear,2);
EndNear= 0;

block_partition(1,1)

Tfar = Tfar(1:EndFar,:);
Tnear = Tnear(1:EndNear,:);

    function block_partition(i,j)
        sigma = ReceiverTree(i);
        tau = SourceTree(j);
        if is_admissible(sigma, tau)
            EndFar = EndFar+1;
            if EndFar > CapacityFar
                CapacityFar = 2*CapacityFar;
                Tfar(CapacityFar,2) = 0;
            end
            Tfar(EndFar,:) = [i, j];
        elseif isempty(sigma.children) && isempty(tau.children)
            EndNear = EndNear+1;
            if EndNear > CapacityNear
                CapacityNear = 2*CapacityNear;
                Tnear(CapacityNear,2) = 0;
            end
            Tnear(EndNear,:) = [i, j];
        else
            Ds = norm(diff(sigma.bb));
            Dt = norm(diff(tau.bb));
            if (Ds > Dt && ~isempty(sigma.children)) || isempty(tau.children) % divide sigma
                s = sigma.children;
                for k = 1 : length(s)
                    block_partition(s(k),j);
                end
            else % divide tau
                t = tau.children;
                for k = 1 : length(t)
                    block_partition(i,t(k));
                end
            end
        end
    end
end

function b = is_admissible(C1, C2, eta)
if nargin < 3
    eta = .8;
end
bb1 = C1.bb;
bb2 = C2.bb;
d1 = norm(diff(bb1));
d2 = norm(diff(bb2));
dist = max(norm(mean(bb1)-mean(bb2))-(d1+d2)/2, 0);
b = min(d1, d2) < eta * dist;
end
