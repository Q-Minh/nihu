classdef Mesh < handle
    properties (SetAccess = private)
        Tolerance
        PointSet
        ElementSets
    end
    
    methods
        function obj = Mesh(coordSys)
            if nargin < 1
                coordSys = DefCoordinateSystem.C3;
            end
            
            if ~isa(coordSys, 'CoordinateSystem')
                error('NiHu:PointSet:invalid_argument', ...
                    'argument ''coordSys'' must be a ''CoordinateSystem'' instance');
            end
            
            obj.PointSet = PointSet(coordSys);
            obj.ElementSets = containers.Map('keyType', 'char', 'valueType', 'any');
            obj.Tolerance = 1e-10;
        end
        
        function addElementSet(obj, elemSet, id)
            if nargin < 3
                id = obj.generateId();
            end
            
            if ~isa(elemSet, 'ElementSet')
                error('NiHu:PointSet:invalid_argument', ...
                    'argument ''elemSet'' must be an ''ElementSet'' instance');
            end
            
            newPointIds = obj.PointSet.merge(elemSet.PointSet, obj.Tolerance);
            oldPointIds = elemSet.PointSet.PointIds;
            elemIds = elemSet.ElemIds;
            elemMatrix = elemSet.ElemMatrix;

            em = elemMatrix(:,2:end);
            [a, b] = ismember(em(:), oldPointIds);
            em(a) = newPointIds(b(a));
            elemMatrix(:,2:end) = em;
            
            obj.ElementSets(id) = ElementSet(obj.PointSet, elemIds, elemMatrix);
        end
    end
    
    methods (Access = private)
        function id = generateId(obj)
            id = sprintf('Set%d', obj.ElementSets.Count()+1);
        end
    end
end
