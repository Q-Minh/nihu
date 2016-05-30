classdef ShapeSet < handle
    %SHAPESET Shape function set representations
    
    % Last modified: 2015.05.10. FP. Introduced Bending shape sets
    
    properties (SetAccess = immutable)
        Id      % unique identifier of the shape function set
        
        Domain  % the shape set's intrinsic domain
        Nodes   % the nodal locations of the shape set
        
        ConstructionMethod
        Args
    end % of immutable properties
    
    methods (Access = private)
        function res = compute_shape_functions(obj)
            fprintf(1, 'Constructing shape functions for id %d\n', obj.Id);
            
            switch obj.ConstructionMethod
                case 'lagrange'
                    base = obj.Args{1};
                    [res.N, g] = lagrangePoly(base, obj.Nodes);
                case 'lagrange-infinite'
                    base = obj.Args{1};
                    [res.N, g] = infiniteLagrangePoly(base, obj.Nodes);
                case 'bending'
                    [res.N, g] = bendingPoly(obj.Domain);
                otherwise
                    error('NiHu:ShapeSet:Arg',...
                        'unknown construction method %s', obj.ConstructionMethod);
            end
            
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
        function obj = ShapeSet(id, domain, nodes, type, varargin)
            obj.Id = id;
            obj.Domain = domain;
            obj.Nodes = nodes;
            obj.ConstructionMethod = type;
            obj.Args = varargin;
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
            
            if numel(ComputedShapeFunctions) < obj.Id ||...
                    isempty(ComputedShapeFunctions{obj.Id})
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
            obj(id == 010, 1) = ShapeSet.ConstantPoint;
            obj(id == 120, 1) = ShapeSet.ConstantLine;
            obj(id == 121, 1) = ShapeSet.LinearLine;
            obj(id == 122, 1) = ShapeSet.QuadraticLine;
            obj(id == 1010, 1) = ShapeSet.InfiniteLine;
            obj(id == 230, 1) = ShapeSet.ConstantTria;
            obj(id == 231, 1) = ShapeSet.LinearTria;
            obj(id == 232, 1) = ShapeSet.QuadraticTria;
            obj(id == 240, 1) = ShapeSet.ConstantQuad;
            obj(id == 241, 1) = ShapeSet.LinearQuad;
            obj(id == 1121, 1) = ShapeSet.InfiniteLinearQuad;
            obj(id == 242, 1) = ShapeSet.QuadraticQuad;
            obj(id == 340, 1) = ShapeSet.ConstantTetra;
            obj(id == 341, 1) = ShapeSet.LinearTetra;
            obj(id == 360, 1) = ShapeSet.ConstantPenta;
            obj(id == 361, 1) = ShapeSet.LinearPenta;
            obj(id == 1231, 1) = ShapeSet.InfiniteLinearPenta;
            obj(id == 380, 1) = ShapeSet.ConstantHexa;
            obj(id == 381, 1) = ShapeSet.LinearHexa;
            obj(id == 1241, 1) = ShapeSet.InfiniteLinearHexa;
            obj(id == 2223, 1) = ShapeSet.BendingLine;
            obj(id == 2333, 1) = ShapeSet.BendingTria;
            obj(id == 2343, 1) = ShapeSet.BendingQuad;
        end
        
        function cid = constant(id)
            data = [
                120 120
                121 120
                122 120
                230 230
                231 230
                232 230
                240 240
                241 240
                242 240
                340 340
                341 340
                342 340
                360 360
                361 360
                380 380
                381 380
                ];
            [a, n] = ismember(id,data(:,1));
            if any(~a)
                error('NiHu:ShapeSetconstant',...
                    'Some shape set IDs could not be transformed to constant');
            end
            cid = data(n,2);
        end
    end % of static methods
    
    enumeration
        ConstantPoint(010, Domain.Point,...
            Domain.Point.Center,...
            'lagrange', [])
        ConstantLine(120, Domain.Line,...
            Domain.Line.Center,...
            'lagrange', [0]) %#ok<NBRAK>
        LinearLine(121, Domain.Line,...
            Domain.Line.CornerNodes,...
            'lagrange', [0; 1])
        QuadraticLine(122, Domain.Line,...
            [-1; 0; 1],...
            'lagrange', [0; 1; 2])
        InfiniteLine(1010, Domain.Line,...
            [-1; 0], 'lagrange-infinite', [0]); %#ok<NBRAK>
        ConstantTria(230, Domain.Tria,...
            Domain.Tria.Center,...
            'lagrange', [0 0])
        LinearTria(231, Domain.Tria,...
            Domain.Tria.CornerNodes,...
            'lagrange', [0 0; 1 0; 0 1])
        QuadraticTria(232, Domain.Tria,...
            [0 0; .5 0; 1 0; .5 .5; 0 1; 0 .5],...
            'lagrange', [0 0; 1 0; 0 1; 2 0; 1 1; 0 2])
        ConstantQuad(240, Domain.Quad,...
            Domain.Quad.Center,...
            'lagrange', [0 0])
        LinearQuad(241, Domain.Quad,...
            Domain.Quad.CornerNodes,...
            'lagrange', [0 0; 1 0; 0 1; 1 1])
        QuadraticQuad(242, Domain.Quad,...
            [-1 -1; 0 -1; 1 -1; 1 0; 1 1; 0 1; -1 1; -1 0],...
            'lagrange', [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 2 1; 1 2])
        InfiniteLinearQuad(1121, Domain.Quad,...
            [-1 -1; 1 -1; 1 0; -1 0], 'lagrange-infinite', [0; 1]);
        ConstantTetra(360, Domain.Tetra,...
            Domain.Tetra.Center,...
            'lagrange', [0 0 0])
        LinearTetra(341, Domain.Tetra,...
            Domain.Tetra.CornerNodes, ...
            'lagrange', [0 0 0; 1 0 0; 0 1 0; 0 0 1]);
        ConstantPenta(360, Domain.Penta,...
            Domain.Penta.Center, ...
            'lagrange', [0 0 0]);
        LinearPenta(361, Domain.Penta,...
            Domain.Penta.CornerNodes, ...
            'lagrange', [0 0 0; 1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1]);
        InfiniteLinearPenta(1231, Domain.Penta,...
            [0 0 -1; 1 0 -1; 0 1 -1; 0 0 0; 1 0 0; 0 1 0],...
            'lagrange-infinite', [0 0; 1 0; 0 1]);
        ConstantHexa(380, Domain.Hexa,...
            Domain.Hexa.Center,...
            'lagrange', [0 0 0])
        LinearHexa(381, Domain.Hexa,...
            Domain.Hexa.CornerNodes,...
            'lagrange', [
            0 0 0; 1 0 0; 0 1 0; 0 0 1
            1 1 0; 0 1 1; 1 0 1; 1 1 1
            ]);
        InfiniteLinearHexa(1241, Domain.Hexa,...
            [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1;
            -1 -1 0; 1 -1 0; 1 1 0; -1 1 0], 'lagrange-infinite',...
            [0 0; 1 0; 0 1; 1 1]);
        BendingLine(2223, Domain.Quad,Domain.Quad.CornerNodes,...
            'bending');
        BendingQuad(2343, Domain.Quad,Domain.Quad.CornerNodes,...
            'bending');
        BendingTria(2333, Domain.Tria,Domain.Tria.CornerNodes,...
            'bending');
    end % of enumeration
end % of class
