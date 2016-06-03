classdef Mesh < handle
    properties (SetAccess = private)
        Tolerance
        PointSet
        ElementSets
        CoordinateSystem
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
            
            obj.CoordinateSystem = coordSys;
            obj.PointSet = PointSet(obj.CoordinateSystem);
            obj.ElementSets = containers.Map('keyType', 'char', 'valueType', 'any');
            obj.Tolerance = 1e-10;
        end % of Constructor
        
        function addElementSet(obj, elemSet, id)
            if nargin < 3
                id = obj.generateElementSetId();
            end
            
            if ~isa(elemSet, 'ElementSet')
                error('NiHu:PointSet:invalid_argument', ...
                    'argument ''elemSet'' must be an ''ElementSet'' instance');
            end
            
            newPointIds = obj.PointSet.merge(elemSet.PointSet);
            oldPointIds = elemSet.PointSet.PointIds;
            elemIds = elemSet.ElemIds;
            elemMatrix = elemSet.ElemMatrix;

            em = elemMatrix(:,2:end);
            [a, b] = ismember(em(:), oldPointIds);
            em(a) = newPointIds(b(a));
            elemMatrix(:,2:end) = em;
            
            obj.ElementSets(id) = ElementSet(obj.PointSet, elemIds, elemMatrix);
        end
        
        function translate(obj, v)
            obj.PointSet.translate(v);
        end % of function translate
    end
    
    methods (Access = private)
        function id = generateElementSetId(obj)
            id = sprintf('Set%d', obj.ElementSets.Count()+1);
        end
    end
end
