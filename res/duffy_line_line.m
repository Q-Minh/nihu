function [Xi, W] = duffy_line_line(N, type)

[x, w] = gaussquad(N, 0, 1);
[x1, x2] = ndgrid(x, x);
w = w(:) * w(:).';
x = [x1(:) x2(:)];
w = w(:);

switch type
    case 'edge'
        Xi = [
            x(:,1), x(:,1) .* x(:,2)
            x(:,1) .* x(:,2), x(:,2)
            ];
        W = [
            w .* x(:,1)
            w .* x(:,2)
            ];
end

end