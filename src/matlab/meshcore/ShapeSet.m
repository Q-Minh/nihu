classdef ShapeSet < handle
    %SHAPESET Shape function set representations
    
    % Last modified: 2015.05.10. FP. Introduced Bending shape sets
    
    properties (SetAccess = immutable)
        Type
        Id      % unique identifier of the shape function set
        
        Domain  % the shape set's intrinsic domain
        Nodes   % the nodal locations of the shape set
        
        Args
    end % of immutable properties
    
    methods (Access = private)
        
        function res = compute_shape_functions(obj)
            fprintf(1, 'Constructing shape functions for id %d\n', obj.Id);
            
            [res.N, g] = obj.Type.construct(obj.Nodes);
            
            nSh = length(res.N);
            d = obj.Domain.Space.Dimension;
            
            res.dN = sym('dN', [d, nSh]);
            res.ddN = sym('ddN', [d*d, nSh]);
            for j = 1 : d
                res.dN(j,:) = simplify(diff(res.N, g(j)));
                for k = 1 : d
                    res.ddN((j-1)*d+k,:) = simplify(diff(res.dN(j,:), g(k)));
                end
            end
            
            res.NFunc = matlabFunction(res.N, 'vars', {g});
            res.dNFunc = cell(nSh, d);
            res.ddNFunc = cell(nSh, d, d);
            for i = 1 : nSh
                for j = 1 : d
                    res.dNFunc{i,j} = matlabFunction(res.dN(j,i), 'vars', {g});
                    for k = 1 : d
                        res.ddNFunc{i,j,k} = matlabFunction(simplify(res.ddN((j-1)*d+k,i)), 'vars', {g});
                    end
                end
            end
        end % of function compute_shape_functions
    end % of private methods
    
    methods (Access = public)
        function obj = ShapeSet(domain, nodes, type)
            obj.Type = type;
            obj.Id = uint32(type.Id * 100 + domain.Id);
            obj.Domain = domain;
            obj.Nodes = nodes;
        end % of constructor
        
        function [N, dN, ddN] = eval(obj, xi)
            %EVAL evaluate shape functions
            % [N, dN, ddN] = EVAL(obj, xi) evaluates the shape functions in
            % the locations XI
            
            % compute symbolic shape function expressions
            persistent ComputedShapeFunctions
            
            if isempty(ComputedShapeFunctions)
                try
                    path = fileparts(mfilename('fullpath'));
                    load(fullfile(path, 'ComputedShapeFunctions.mat'));
                catch
                    ComputedShapeFunctions = containers.Map('keyType', 'uint32', 'valueType', 'any');
                end
            end
            
            if ~ComputedShapeFunctions.isKey(obj.Id)
                shfun = obj.compute_shape_functions();
                ComputedShapeFunctions(obj.Id) = shfun;
                path = fileparts(mfilename('fullpath'));
                save(fullfile(path, 'ComputedShapeFunctions.mat'), 'ComputedShapeFunctions');
            else
                shfun = ComputedShapeFunctions(obj.Id);
            end
            
            q = size(xi,1); % number of locations
            d = size(xi,2); % number of dimensions
            nSh = length(shfun.N);  % number of shape functions
            
            N = shfun.NFunc(xi);
            if size(N,1) == 1 && q > 1
                N = repmat(N, q, 1);
            end
            
            % evaluate shape function derivatives
            if nargout == 1
                return
            end
            
            dN = zeros(q, nSh, d);
            for i = 1 : nSh
                for k = 1 : d
                    dN(:,i,k) = shfun.dNFunc{i,k}(xi);
                end
            end
            
            % evaluate shape function second derivatives
            if nargout == 2
                return
            end
            
            ddN = zeros(q, nSh, d, d);
            for i = 1 : nSh
                for j = 1 : d
                    for k = 1 : d
                        ddN(:,i,j,k) = shfun.ddNFunc{i,j,k}(xi);
                    end
                end
            end
        end % of function eval
    end % of public methods
    
    methods (Static = true)
        function obj = fromId(id)
            obj = ShapeSet.empty(size(id,1),0);
            obj(id == ShapeSet.ConstantPoint.Id, 1) = ShapeSet.ConstantPoint;
            obj(id == ShapeSet.ConstantLine.Id, 1) = ShapeSet.ConstantLine;
            obj(id == ShapeSet.LinearLine.Id, 1) = ShapeSet.LinearLine;
            obj(id == ShapeSet.QuadraticLine.Id, 1) = ShapeSet.QuadraticLine;
            obj(id == ShapeSet.InfiniteLine.Id, 1) = ShapeSet.InfiniteLine;
            obj(id == ShapeSet.ConstantTria.Id, 1) = ShapeSet.ConstantTria;
            obj(id == ShapeSet.LinearTria.Id, 1) = ShapeSet.LinearTria;
            obj(id == ShapeSet.QuadraticTria.Id, 1) = ShapeSet.QuadraticTria;
            obj(id == ShapeSet.ConstantQuad.Id, 1) = ShapeSet.ConstantQuad;
            obj(id == ShapeSet.LinearQuad.Id, 1) = ShapeSet.LinearQuad;
            obj(id == ShapeSet.InfiniteLinearQuad.Id, 1) = ShapeSet.InfiniteLinearQuad;
            obj(id == ShapeSet.QuadraticQuad.Id, 1) = ShapeSet.QuadraticQuad;
            obj(id == ShapeSet.QuadraticQuadMid.Id, 1) = ShapeSet.QuadraticQuadMid;
            obj(id == ShapeSet.ConstantTetra.Id, 1) = ShapeSet.ConstantTetra;
            obj(id == ShapeSet.LinearTetra.Id, 1) = ShapeSet.LinearTetra;
            obj(id == ShapeSet.ConstantPenta.Id, 1) = ShapeSet.ConstantPenta;
            obj(id == ShapeSet.LinearPenta.Id, 1) = ShapeSet.LinearPenta;
            obj(id == ShapeSet.InfiniteLinearPenta.Id, 1) = ShapeSet.InfiniteLinearPenta;
            obj(id == ShapeSet.ConstantHexa.Id, 1) = ShapeSet.ConstantHexa;
            obj(id == ShapeSet.LinearHexa.Id, 1) = ShapeSet.LinearHexa;
            obj(id == ShapeSet.InfiniteLinearHexa.Id, 1) = ShapeSet.InfiniteLinearHexa;
            obj(id == ShapeSet.BendingLine.Id, 1) = ShapeSet.BendingLine;
            obj(id == ShapeSet.BendingTria.Id, 1) = ShapeSet.BendingTria;
            obj(id == ShapeSet.BendingQuad.Id, 1) = ShapeSet.BendingQuad;
        end
        
        function nsetId = constant(domainId)
            %CONSTANT assign constant shape set to a domain
            data = [
                Domain.Line.Id ShapeSet.ConstantLine.Id
                Domain.Tria.Id ShapeSet.ConstantTria.Id
                Domain.Quad.Id ShapeSet.ConstantQuad.Id
                Domain.Tetra.Id ShapeSet.ConstantTetra.Id
                Domain.Penta.Id ShapeSet.ConstantPenta.Id
                Domain.Hexa.Id ShapeSet.ConstantHexa.Id
                ];
            [a, n] = ismember(domainId, data(:,1));
            if any(~a)
                error('NiHu:ShapeSetconstant',...
                    'Some shape set IDs could not be transformed to constant');
            end
            nsetId = data(n,2);
        end
    end % of static methods
    
    enumeration
        ConstantPoint(Domain.Point,...
            Domain.Point.Center,...
            LagrangeShapeSet(0, []))
        ConstantLine(Domain.Line,...
            Domain.Line.Center,...
            LagrangeShapeSet(0, [0])) %#ok<NBRAK>
        LinearLine(Domain.Line,...
            Domain.Line.CornerNodes,...
            LagrangeShapeSet(1, [0; 1]))
        QuadraticLine(Domain.Line,...
            [-1; 0; 1],...
            LagrangeShapeSet(2, [0; 1; 2]))
        InfiniteLine(Domain.Line,...
            [-1; 0], InfiniteLagrangeShapeSet(0, [0])); %#ok<NBRAK>
        ConstantTria(Domain.Tria,...
            Domain.Tria.Center,...
            LagrangeShapeSet(0, [0 0]))
        LinearTria(Domain.Tria,...
            Domain.Tria.CornerNodes,...
            LagrangeShapeSet(1, [0 0; 1 0; 0 1]))
        QuadraticTria(Domain.Tria,...
            [0 0; .5 0; 1 0; .5 .5; 0 1; 0 .5],...
            LagrangeShapeSet(2, [0 0; 1 0; 0 1; 2 0; 1 1; 0 2]))
        ConstantQuad(Domain.Quad,...
            Domain.Quad.Center,...
            LagrangeShapeSet(0, [0 0]))
        LinearQuad(Domain.Quad,...
            Domain.Quad.CornerNodes,...
            LagrangeShapeSet(1, [0 0; 1 0; 0 1; 1 1]))
        QuadraticQuad(Domain.Quad,...
            [-1 -1; 0 -1; 1 -1; 1 0; 1 1; 0 1; -1 1; -1 0],...
            LagrangeShapeSet(2, [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 2 1; 1 2]))
        QuadraticQuadMid(Domain.Quad,...
            [-1 -1; 0 -1; 1 -1; 1 0; 1 1; 0 1; -1 1; -1 0; 0 0],...
            LagrangeShapeSet(2, [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 2 1; 1 2; 2 2]))
        InfiniteLinearQuad(Domain.Quad,...
            [-1 -1; 1 -1; 1 0; -1 0], InfiniteLagrangeShapeSet(1, [0; 1]));
        ConstantTetra(Domain.Tetra,...
            Domain.Tetra.Center,...
            LagrangeShapeSet(0, [0 0 0]))
        LinearTetra(Domain.Tetra,...
            Domain.Tetra.CornerNodes, ...
            LagrangeShapeSet(1, [0 0 0; 1 0 0; 0 1 0; 0 0 1]));
        ConstantPenta(Domain.Penta,...
            Domain.Penta.Center, ...
            LagrangeShapeSet(0, [0 0 0]));
        LinearPenta(Domain.Penta,...
            Domain.Penta.CornerNodes, ...
            LagrangeShapeSet(1, [0 0 0; 1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1]));
        InfiniteLinearPenta(Domain.Penta,...
            [0 0 -1; 1 0 -1; 0 1 -1; 0 0 0; 1 0 0; 0 1 0],...
            InfiniteLagrangeShapeSet(1, [0 0; 1 0; 0 1]));
        ConstantHexa(Domain.Hexa,...
            Domain.Hexa.Center,...
            LagrangeShapeSet(0, [0 0 0]))
        LinearHexa(Domain.Hexa,...
            Domain.Hexa.CornerNodes,...
            LagrangeShapeSet(1, [
            0 0 0; 1 0 0; 0 1 0; 0 0 1
            1 1 0; 0 1 1; 1 0 1; 1 1 1
            ]));
        InfiniteLinearHexa(Domain.Hexa,...
            [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1;
            -1 -1 0; 1 -1 0; 1 1 0; -1 1 0], InfiniteLagrangeShapeSet(1,...
            [0 0; 1 0; 0 1; 1 1]));
        BendingLine(Domain.Line, Domain.Line.CornerNodes,...
            BendingShapeSet(3, Domain.Line));
        BendingTria(Domain.Tria, Domain.Tria.CornerNodes,...
            BendingShapeSet(3, Domain.Tria));
        BendingQuad(Domain.Quad, Domain.Quad.CornerNodes,...
            BendingShapeSet(3, Domain.Quad));
    end % of enumeration
end % of class
