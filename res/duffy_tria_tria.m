function [Xi, W] = duffy_tria_tria(N, type)

% generate raw cube quadrature
[x, w] = gaussquad(N, 0, 1);
[x1,x2,x3,x4] = ndgrid(x,x,x,x);
x = [x1(:) x2(:) x3(:) x4(:)];
w = w * w';
w = w(:);
w = w * w';
w = w(:);

% transform into quadrature on triangle simplex
Xi = [];
W = [];
switch lower(type)
    case 'face'
        for i = 1 : 3
            [xi, ww] = trans_tria_face(i, x, w);
            Xi = [Xi; xi];
            W = [W; ww];
        end
    case 'edge'
        for i = 1 : 6
            [xi, ww] = trans_tria_edge(i, x, w);
            Xi = [Xi; xi];
            W = [W; ww];
        end
    case 'corner'
        [Xi, W] = trans_tria_corner(x, w);
end

% Apply transform to our standard triangle
T = [
    1 0 0 0
    -1 1 0 0
    0 0 1 0
    0 0 -1 1
    ];
Xi = Xi * T;

end

function [xi, ww] = trans_tria_face(n, x, w)
switch n
    case 1
        mu1 = x(:,1);
        mu2 = x(:,1) .* x(:,2);
        xi1 = (1-mu1) .* x(:,3);
        xi2 = xi1 .* x(:,4);
        J = (1-mu1) .* xi1;
    case 2
        mu1 = x(:,1) .* x(:,2);
        mu2 = x(:,1) .* (x(:,2) - 1);
        xi1 = (1-mu1+mu2) .* x(:,3) - mu2;
        xi2 = (xi1+mu2) .* x(:,4) - mu2;
        J = (1-mu1+mu2) .* (xi1+mu2);
    case 3
        mu1 = x(:,1) .* x(:,2);
        mu2 = x(:,1);
        xi1 = (1-mu2) .* x(:,3) + mu2 - mu1;
        xi2 = (xi1-mu2+mu1) .* x(:,4);
        J = (1-mu2) .* (xi1-mu2+mu1);
end

J = J .* x(:,1);

xi = [
    xi1 xi2 xi1+mu1 xi2+mu2
    xi1+mu1 xi2+mu2 xi1 xi2
    ];

ww = [w .* J; w .* J];

end



function [xi, ww] = trans_tria_edge(n, x, w)
switch n
    case 1
        mu1 = -x(:,1) .* x(:,2);
        mu2 = -x(:,1) .* x(:,2) .* x(:,3);
        xi1 = (1-x(:,1)) .* x(:,4) + x(:,1);
        xi2 = x(:,1) .* (1-x(:,2)+x(:,2).*x(:,3));
    case 2
        mu1 = x(:,1) .* x(:,2);
        mu2 = x(:,1) .* x(:,2) .* x(:,3);
        xi1 = (1-x(:,1)) .* x(:,4) + x(:,1) .* (1-x(:,2));
        xi2 = x(:,1) .* (1-x(:,2));
    case 3
        mu1 = -x(:,1) .* x(:,2) .* x(:,3);
        mu2 = x(:,1) .* x(:,2) .* (1 - x(:,3));
        xi1 = (1-x(:,1)) .* x(:,4) + x(:,1);
        xi2 = x(:,1) .* (1-x(:,2));
    case 4
        mu1 = x(:,1) .* x(:,2) .* x(:,3);
        mu2 = x(:,1) .* x(:,2) .* (x(:,3) - 1);
        xi1 = (1-x(:,1)) .* x(:,4) + x(:,1) .* (1-x(:,2).*x(:,3));
        xi2 = x(:,1) .* (1-x(:,2).*x(:,3));
    case 5
        mu1 = -x(:,1) .* x(:,2) .* x(:,3);
        mu2 = -x(:,1) .* x(:,2);
        xi1 = (1-x(:,1)) .* x(:,4) + x(:,1);
        xi2 = x(:,1);
    case 6
        mu1 = x(:,1) .* x(:,2) .* x(:,3);
        mu2 = x(:,1) .* x(:,2);
        xi1 = (1-x(:,1)) .* x(:,4) + x(:,1) .* (1-x(:,2).*x(:,3));
        xi2 = x(:,1) .* (1-x(:,2));
end

xi = [xi1 xi2 xi1+mu1 xi2+mu2];
J = x(:,2) .* x(:,1).^2 .* (1-x(:,1));
ww = w .* J;

end


function [xi, ww] = trans_tria_corner(x, w)
xi1 = x(:,1);
xi2 = xi1.*x(:,2);
eta1 = xi1.*x(:,3);
eta2 = eta1.*x(:,4);
J = x(:,3).*x(:,1).^3;

xi = [
    xi1 xi2 eta1 eta2
    eta1 eta2 xi1 xi2
    ];

ww = [w .* J; w .* J];

end
