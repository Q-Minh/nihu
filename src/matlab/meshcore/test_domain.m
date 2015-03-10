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
    fprintf('Testing Domain %s...\n', d.char());

    %% check domain volume with 0-th order numerical quadrature
    [~, w] = gaussian_quadrature(d, 1);
    err = abs(sum(w)/d.Volume-1);
    if (err > 1e-3)
        error('Domain volume error: %s', d.char());
    end
    
    %% check center with numerical integration
    [xi, w] = gaussian_quadrature(d, 2);
    c_num = w.' * xi / sum(w);
    err = norm(c_num-d.Center)/norm(d.Center);
    if (err > 1e-3)
        error('Domain center error: %s', d.char());
    end
    
    %% 
    disp('Tests OK');
end
