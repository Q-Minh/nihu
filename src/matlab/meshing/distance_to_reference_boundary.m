function rho = distance_to_reference_boundary(xi, theta, domain)
switch domain
    case Domain.Tria
        rho(:,1) = -xi(1) ./ cos(theta);
        rho(:,2) = -xi(2) ./ sin(theta);
        rho(:,3) = (1-sum(xi)) ./ (sin(theta) + cos(theta));
    case Domain.Quad
        rho(:,1) = (-1 - xi(1)) ./ sin(theta);
        rho(:,2) = (1 - xi(2)) ./ cos(theta);
        rho(:,3) = (1 - xi(1)) ./ sin(theta);
        rho(:,4) = (-1 - xi(2)) ./ cos(theta);
end
rho(rho < 0) = Inf;
rho = min(rho, [], 2);
end
