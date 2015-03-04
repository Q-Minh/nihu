classdef ShapeSet
    properties
        Domain  % the shape set's intrinsic domain
        Nodes   % the nodal locations of the shape set
        
        N       % symbolic expressions of the shape functions
        dN      % symbolic expressions of the shape functions's derivatives
        
        NFunc   % Matlab function evaluating the shape functions
        dNFunc  % Matlab function evaluating the shape function derivatives
    end
    
    methods
        function obj = ShapeSet(domain, nodes, type, varargin)
            obj.Domain = domain;
            obj.Nodes = nodes;
            
            switch type
                case 'lagrange'
                    base = varargin{1};
                    [obj.N, obj.dN, obj.NFunc, obj.dNFunc] =...
                        lagrangePoly(base, nodes);
            end
        end % of constructor
        
        function [N, dN] = eval(obj, xi)
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
        end % of function eval
    end
    
    enumeration
        LinearLine(Domain.Line,...
            [-1; 1],...
            'lagrange', [0; 1])
        QuadraticLine(Domain.Line,...
            [-1; 0; 1],...
            'lagrange', [0; 1; 2])
        LinearTria(Domain.Tria,...
            [0 0; 1 0; 0 1],...
            'lagrange', [0 0; 1 0; 0 1])
        QuadraticTria(Domain.Tria,...
            [0 0; .5 0; 1 0; .5 .5; 0 1; 0 .5],...
            'lagrange', [0 0; 1 0; 0 1; 2 0; 1 1; 0 2])
        LinearQuad(Domain.Quad,...
            [-1 -1; 1 -1; 1 1; -1 1],...
            'lagrange', [0 0; 1 0; 0 1; 1 1])
        QuadraticQuad(Domain.Quad,...
            [-1 -1; 0 -1; 1 -1; 1 0; 1 1; 0 1; -1 1; -1 0],...
            'lagrange', [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 2 1; 1 2])
    end
end
