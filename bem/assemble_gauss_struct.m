function g = assemble_gauss_struct(nvert, nnod)
%ASSEMBLE_GAUSS_STRUCT
%  G = ASSEMBLE_GAUSS_STRUCT Assembles the Gaussian quadrature structure
%  used by the functions BEMHG.
%  This is an array of structures with fields
%    w    : 2D Gaussian weights
%    N    : Linear shape functions
%    Nxi  : Derivative of shape functions in the first direction
%    Neta : Derivative of shape functions in the second direction
%  The array G is assembled as follows: The last entries describe the
%  quadrature used for singular elements. The 1., 2., 3.  values are the near
%  field, mid field and far field quadrature structures.
%
% See also: gaussquad2, shapefun

% Peter Fiala
% Last modified: 14.11.2009.

nNod = length(nnod);
c = cell(nNod,1);
g = struct('w', c, 'N', c, 'Nxi', c, 'Neta', c, 'xi', c);
for iNod = 1 : nNod
    [g(iNod).xi, g(iNod).w] = gaussquad2(nnod(iNod),nvert);
    [g(iNod).N, dN] = shapefun(g(iNod).xi, 20+nvert);
    g(iNod).Nxi = dN(:,:,1);
    g(iNod).Neta = dN(:,:,2);
end
end
