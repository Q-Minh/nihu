classdef PolynomialShapeSet < GenericShapeSet
    methods (Abstract = true)
        [N, g] = construct(base, nodes);
    end
    
    properties (SetAccess = immutable)
        Order
        Base
        Id
    end
    
    methods
        function obj = PolynomialShapeSet(order, base)
            obj.Order = order;
            obj.Base = base;
            obj.Id = 10 * obj.FamilyId + obj.Order;
        end
    end % of methods
end % of PolynomialShapeSet class
