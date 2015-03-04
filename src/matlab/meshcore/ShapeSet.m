classdef ShapeSet
    properties
        Domain
        Nodes
        
        NFunc;
        dNFunc;
    end
    
    methods
        function obj = ShapeSet(domain, nodes, type, varargin)
            obj.Domain = domain;
            obj.Nodes = nodes;
            
            switch type
                case 'lagrange'
                    [~, obj.NFunc, obj.dNFunc] = lagrangePoly(varargin{1}, nodes);
            end
        end
        
        function [N, dN] = eval(obj, xi)
            N = obj.NFunc(xi);
            if size(N,1) == 1 && size(xi,1) > 1
                N = repmat(N, size(xi,1), 1);
            end
            
            dN = zeros(size(xi,1), size(obj.Nodes,1), size(xi,2));
            for k = 1 : size(xi,2)
                dn = obj.dNFunc{k}(xi);
                if size(dn,1) == 1 && size(xi,1) > 1
                    dN(:,:,k) = repmat(dn, size(xi,1), 1);
                else
                    dN(:,:,k) = dn;
                end
            end
        end
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
        LinearQuad(Domain.Quad,...
            [-1 -1; 1 -1; 1 1; -1 1],...
            'lagrange', [0 0; 1 0; 0 1; 1 1])
    end
end
