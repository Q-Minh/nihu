clear;
clc;

domains = [
    Domain.Line
    Domain.Tria
    Domain.Quad
    Domain.Tetra
    Domain.Penta
    Domain.Hexa
];

for i = 1 : length(domains)
    d = domains(i);
    disp(d.char());
    [xi, w] = gaussian_quadrature(d, 4);
    disp(sum(w));
    disp(w.' * xi / sum(w));
end
