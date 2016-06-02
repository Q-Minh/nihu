classdef ElementSet < handle
    properties (SetAccess = private)
        ElementSetId
        PointSet
        ElemIds
        ElemMatrix % ElemId, ShapeSetId, node Ids ... 0 padding
    end
    
    methods
        function obj = ElementSet(id, pointSet)
            % error check
            obj.ElementSetId = uint32(id);
            obj.PointSet = pointSet;
            obj.ElemIds = uint32(zeros(0,1));
            obj.ElemMatrix = uint32(zeros(0,1));
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
