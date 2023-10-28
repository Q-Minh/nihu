function [h, hn, hi, ha] = plot_cluster(CLevel, cluster)

% plot cluster
h = plot_box(CLevel.coord(cluster,:), CLevel.diameter, 'g');

% plot near field
[~,~,near] = find(CLevel.nearfield(cluster,:));
near = setdiff(near, cluster);
hn = plot_box(CLevel.coord(near,:), CLevel.diameter, 'r');

% plot interaction list
[~,~,il] = find(CLevel.interlist(cluster,:));
il = setdiff(il, cluster);
hi = plot_box(CLevel.coord(il,:), CLevel.diameter, 'b');

% plot everybody
al = 1 : size(CLevel.coord,1);
al = setdiff(al, il);
al = setdiff(al, near);
al = setdiff(al, cluster);
ha = plot_box(CLevel.coord(al,:), CLevel.diameter, 'w');

end

function h = plot_box(x0, D, col)

dim = size(x0,2);

if (dim == 3)
    p = [
        1 2 3 4
        5 6 7 8
        1 2 6 5
        2 3 7 6
        3 4 8 7
        4 1 5 8
        ];
    
    q = [
        -1 -1 -1
        1 -1 -1
        1 1 -1
        -1 1 -1
        -1 -1 1
        1 -1 1
        1 1 1
        -1 1 1
        ] * D/2;
    x = q(:,1);
    y = q(:,2);
    z = q(:,3);
    h = patch(x(p')+x0(1), y(p')+x0(2), z(p')+x0(3), ones(size(p')));
    
elseif dim == 2
    q = [
        -1 -1
        1 -1
        1 1
        -1 1
        ] * D/2;

    x = q(:,1);
    y = q(:,2);

    h = patch(bsxfun(@plus, x, x0(:,1)'),...
        bsxfun(@plus, y, x0(:,2)'),...
        col);
end

end
