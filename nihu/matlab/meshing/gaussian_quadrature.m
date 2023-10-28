function [xi, w] = gaussian_quadrature(domain, order)
switch domain
    case Domain.Line
        [xi, w] = gaussquad1(order);
    case {Domain.Tria, Domain.Quad}
        nvert = size(domain.CornerNodes,1);
        [xi, w] = gaussquad2(order, nvert);
    case {Domain.Tria, Domain.Quad}
        nvert = size(domain.CornerNodes,1);
        [xi, w] = gaussquad2(order, nvert);
    case {Domain.Tetra, Domain.Penta, Domain.Hexa}
        nvert = size(domain.CornerNodes,1);
        [xi, w] = gaussquad3(order, nvert);
end
end
