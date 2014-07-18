function L = lagpolys(z)
%LAGPOLYS Creates Lagrange polynomials for the given zero points.
%   z      : location of zero points
%   L      : the corresponding Lagrange polynomials as a cell array.

% Peter Rucz
% 2009 November

L = cell(length(z),1);
for c1 = 1:length(z)
    L{c1} = 1;
end
for c1 = 1:length(z);
    ind = [1:c1-1,c1+1:length(z)];
    for c2=1:length(ind)
        L{ind(c2)} = conv(L{ind(c2)},[1 -z(c1)])/prod(z(ind)-z(c1));
    end
end
end