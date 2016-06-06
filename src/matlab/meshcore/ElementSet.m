classdef ElementSet < handle
    properties (SetAccess = private)
        PointSet
        ShapeSet
        ElemIds
        ElemMatrix % node Ids ... 
    end
    
    methods (Static = true)
        function obj = createLine()
            
        end
    end % of Static methods
    
    methods
        function obj = ElementSet(pointSet, shapeSet, elemIds, elemMatrix)
            if nargin < 4
                elemMatrix = uint32(zeros(0,1));
            end

            if nargin < 3
                elemIds = uint32(zeros(0,1));
            end

            % error check
            obj.PointSet = pointSet;
            obj.ShapeSet = shapeSet;
            obj.ElemIds = elemIds;
            obj.ElemMatrix = elemMatrix;
        end
        
        function newId = addElem(obj, elemMatrix, ids)
            if nargin < 3 || isempty(ids)
                n = size(elemMatrix,1);
                newId = obj.generateIds(n);
            else
                if numel(ids) ~= size(elemMatrix,1)
                    error('NiHu:ElementnSet:invalid_argument', ...
                        'Coincident Ids');
                end
                if numel(unique(ids)) ~= numel(ids)
                    error('NiHu:ElementSet:invalid_argument', ...
                        'Id arguments are not unique');
                end
                if any(ismember(ids, obj.ElemIds))
                    error('NiHu:ElementSet:invalid_argument', ...
                        'Coincident Ids');
                end
                newId = ids(:);
            end
            obj.ElemIds = [obj.ElemIds; newId];
            obj.ElemMatrix(end+(1:size(elemMatrix,1)), 1:size(elemMatrix,2)) = uint32(elemMatrix);
        end
        
        function [elements, ids] = getElementsByIndex(obj, idx)
            elements = obj.ElemMatrix(idx,:);
            if nargout > 1
                ids = obj.ElemIds(idx);
            end
        end
        
        function [elements, idx] = getElementsById(obj, ids)
            idx = obj.id2idx(ids);
            elements = getElementsByIndex(obj, idx);
        end
        
        function divide(obj, N)
            [xi, N] = obj.ShapeSet.Domain.divide(N);
            
            xi = reshape(xi(N(:),:), size(N,1), size(N,2), []);
            
            Nc = obj.ShapeSet.eval(obj.ShapeSet.CornerNodes);
            
            for d = 1 : obj.ShapeSet.Domain.Dimension
                x = xi(:, :, d);
                xi_xi(d, :, :) = Nc * x.';
            end
            
            coords = obj.map(xi_xi);
            
        end
        
        function coords = map(obj, xi)
            % TODO: check dimensions
            
            N = obj.ShapeSet.eval(xi);
            
            xyz = obj.PointSet.getCoordinatesById(obj.ElemMatrix(:));
            
            xyz = reshape(xyz, size(obj.ElemMatrix, 1), size(obj.ElemMatrix, 2), []);
            
            coords = zeros(size(obj.PointSet.Coordinates, 2), size(xi, 1), size(obj.ElemMatrix, 1));
            
            for d = 1 : size(obj.PointSet.Coordinates, 2)
                x = xyz(:, :, d);
                coords(d, :, :) = N * x.';
            end
            
        end
        
        
        function translate(obj, v)
            obj.PointSet.translate(v);
        end % of function translate
    end
    
    methods (Access = private)
        function ids = generateIds(obj, n)
            if isempty(obj.ElemIds)
                mx = uint32(0);
            else
                mx = max(obj.ElemIds);
            end
            ids = mx + uint32(1 : n)';
        end
        
        function idx = id2idx(obj, id)
            [ism, idx] = ismember(id, obj.PointIds);
            if ~all(ism)
                error('NiHu:ElementSet:invalid_argument', ...
                    'Id cannot found in ElementSet');
            end
        end
    end
end
