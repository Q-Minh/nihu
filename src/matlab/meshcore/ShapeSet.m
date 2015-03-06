classdef ShapeSet < handle
    properties (SetAccess = immutable)
        Id      % unique identifier of the shape function set
        
        Domain  % the shape set's intrinsic domain
        Nodes   % the nodal locations of the shape set
        
        N       % symbolic expressions of the shape functions
        
        NFunc   % Matlab func evaluating the shape functions
        dNFunc  % Matlab func evaluating the shape function derivatives
        ddNFunc % Matlab func evaluating the shape function 2nd derivatives
    end
    
    methods
        function obj = ShapeSet(id, domain, nodes, type, varargin)
            fprintf(1, 'Constructing shape functions for id %d\n', id);
            
            obj.Id = id;
            obj.Domain = domain;
            obj.Nodes = nodes;
            
            switch type
                case 'lagrange'
                    base = varargin{1};
                    [obj.N, ~, obj.NFunc, obj.dNFunc, obj.ddNFunc] =...
                        lagrangePoly(base, nodes);
            end
        end % of constructor
        
        function [N, dN, ddN] = eval(obj, xi)
            q = size(xi,1);
            d = size(xi,2);
            n = size(obj.Nodes,1);
            
            N = obj.NFunc(xi);
            if size(N,1) == 1 && q > 1
                N = repmat(N, q, 1);
            end
            
            if nargout == 1
                return
            end
            
            dN = zeros(q, n, d);
            for i = 1 : n
                for k = 1 : d
                    dN(:,i,k) = obj.dNFunc{i,k}(xi);
                end
            end
            
            if nargout == 2
                return
            end
            
            ddN = zeros(q, n, d, d);
            for i = 1 : n
                for j = 1 : d
                    for k = 1 : d
                        ddN(:,i,j,k) = obj.ddNFunc{i,j,k}(xi);
                    end
                end
            end
            
        end % of function eval
    end % of public methods
    
    methods (Static = true)
        function obj = fromId(id)
            switch id
                case {12 121}
                    obj = ShapeSet.LinearLine;
                case 122
                    obj = ShapeSet.QuadraticLine;
                case {23 231}
                    obj = ShapeSet.LinearTria;
                case 232
                    obj = ShapeSet.QuadraticTria;
                case {24 241}
                    obj = ShapeSet.LinearQuad;
                case 242
                    obj = ShapeSet.QuadraticQuad;
                case {34 341}
                    obj = ShapeSet.LinearTetra;
                case {36 361}
                    obj = ShapeSet.LinearPenta;
                case {38 381}
                    obj = ShapeSet.LinearHexa;
            end
        end
    end
    
    enumeration
        LinearLine(121, Domain.Line,...
            Domain.Line.CornerNodes,...
            'lagrange', [0; 1])
        QuadraticLine(122, Domain.Line,...
            [-1; 0; 1],...
            'lagrange', [0; 1; 2])
        LinearTria(231, Domain.Tria,...
            Domain.Tria.CornerNodes,...
            'lagrange', [0 0; 1 0; 0 1])
        QuadraticTria(232, Domain.Tria,...
            [0 0; .5 0; 1 0; .5 .5; 0 1; 0 .5],...
            'lagrange', [0 0; 1 0; 0 1; 2 0; 1 1; 0 2])
        LinearQuad(241, Domain.Quad,...
            Domain.Quad.CornerNodes,...
            'lagrange', [0 0; 1 0; 0 1; 1 1])
        QuadraticQuad(242, Domain.Quad,...
            [-1 -1; 0 -1; 1 -1; 1 0; 1 1; 0 1; -1 1; -1 0],...
            'lagrange', [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 2 1; 1 2])
        LinearTetra(341, Domain.Tetra,...
            Domain.Tetra.CornerNodes, ...
            'lagrange', [0 0 0; 1 0 0; 0 1 0; 0 0 1]);
        LinearPenta(361, Domain.Penta,...
            Domain.Penta.CornerNodes, ...
            'lagrange', [0 0 0; 1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1]);
        LinearHexa(381, Domain.Hexa,...
            Domain.Hexa.CornerNodes,...
            'lagrange', [
            0 0 0; 1 0 0; 0 1 0; 0 0 1
            1 1 0; 0 1 1; 1 0 1; 1 1 1
            ]);
    end
end
