function [Xi, W] = DuffyQuad(n, shapeSet, xi0)
%   DUFFYQUAD create Duffy polar transformed quadrature
% [Xi, W] = DUFFYQUAD(N, ID, Xi0) creates an N order Duffy polar transformed
% quadrature over an element of ID, where the singular point is located at
% the local coordinate Xi0

domain = shapeSet.Domain;

order = 2*n-1;
nC = size(domain.CornerNodes, 1);
[xi, w] = gaussquad2(order, nC);
c = domain.CornerNodes;

[N, dN] = shapefun(xi, shapeSet.Id);

Xi = [];
W = [];

C(1,:) = xi0;
C(2,:) = xi0;
for i = 1 : nC
    C(3,:) = c(i,:);
    C(4,:) = c(mod(i+1-1, nC)+1,:);
    
    Xi = [
        Xi
        N * C
        ];
    dx1 = dN(:,:,1) * C;
    dx2 = dN(:,:,2) * C;
    ww = dx1(:,1).*dx2(:,2) - dx2(:,1).*dx1(:,2);
    W = [
        W
        w.*ww
        ];
end

end
