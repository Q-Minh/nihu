function free = unidist(mesh, free)

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

nFree = size(free,1);
fixed = mesh.Nodes(:,2:3);
nFixed = size(fixed,1);

points = [
    free
    fixed
    ];

x = points(:,1);
y = points(:,2);

tri = delaunay(x,y);
edges = [
    tri(:,[1 2])
    tri(:,[2 3])
    tri(:,[1 3])
    ];
edges = unique(sort(edges,2), 'rows');

edg = edges(any(edges <= nFree,2),:);

while true
    x1 = points(edg(:,1),1:2);
    x2 = points(edg(:,2),1:2);
    dvec = x2 - x1;                 % edge distance vectors
    d = sqrt(dot(dvec, dvec, 2));   % edge distances
    d0 = mean(d);                   % mean distance
    gf = 2*dvec.*repmat(1-d0./d,1,2);
    
    f = sum((d-d0).^2);             % cost function
    gradf = zeros(size(points,1),2);
    for iE = 1 : size(edg,1)
        gradf(edg(iE,1),:) = gradf(edg(iE,1),:) - gf(iE,:);
        gradf(edg(iE,2),:) = gradf(edg(iE,2),:) + gf(iE,:);
    end
    
%     mu = d0/2 / max(sqrt(dot(gradf(1:nFree,:),gradf(1:nFree,:),2)));
    mu = .1;
    
    points(1:nFree,:) = points(1:nFree,:) - mu * gradf(1:nFree,:);
    x = points(:,1);
    y = points(:,2);
    
    figure;
    line(x(edges.'), y(edges.'), 'Color', 'red');
    pause;
    
    free = points(1:nFree,:);
end
end
