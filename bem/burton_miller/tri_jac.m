clear;

nodes = [0 0; 1 0; 0 1]*2;

x0 = [1 1]/3;

[xi, w] = gaussquad1(7);


I = 0;

nodepairs = [1 2; 2 3; 3 1];
for ipair = 1:size(nodepairs,1)
    node1 = nodepairs(ipair,1);
    node2 = nodepairs(ipair,2);

    d = nodes(node2,:) - nodes(node1,:);
    L = norm(d);

for ixi = 1:length(xi)
    xip = xi(ixi);
    
    x = nodes(node1,:)*(1-xip)/2 + nodes(node2,:)*(1+xip)/2 - x0;
    r = norm(x);
    sina = sqrt(1 - (dot(x, d) / r / L)^2);
    
    jac = sina/r*L/2;
    I = I+0.5*r*r*w(ixi)*jac;
end
end

