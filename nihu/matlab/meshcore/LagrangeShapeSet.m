classdef LagrangeShapeSet < PolynomialShapeSet
    properties (Constant)
        FamilyId = 0;
    end
    
    methods
        function obj = LagrangeShapeSet(order, base)
            obj = obj@PolynomialShapeSet(order, base);
        end
    end
    
    methods
        function [N, g] = construct(obj, corners)
            
            d = size(corners,2);
            n = size(corners,1);
            
            g = sym('xi', [1 d]);
            s = ones(n,1);
            for j = 1 : d
                s = s .* g(j).^obj.Base(:,j);
            end
            
            q = ones(n,n);
            for k = 1 : n
                for j = 1 : d
                    q(:,k) = q(:,k) .* corners(:,j).^obj.Base(k,j);
                end
            end
            N = simplify(sym(s.' / q));
        end
    end
end
