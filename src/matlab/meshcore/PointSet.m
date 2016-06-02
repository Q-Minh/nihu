classdef PointSet < handle
    properties (SetAccess = private)
        PointSetId
        CoordinateSystem
        PointIds
        Coordinates
    end
    
    methods
        function obj = PointSet(id, coordSys)
            %PointSet constructor
            % obj = PointSet(id, coordSys)
            
            if ~isa(coordSys, 'CoordinateSystem')
                error('NiHu:PointSet:invalid_argument', ...
                    'argument ''coordSys'' must be a ''CoordinateSystem'' instance');
            end
            
            % take Id and CS
            obj.PointSetId = uint32(id);
            obj.CoordinateSystem = coordSys;
            
            % allocate PointIds and Coordinates
            obj.PointIds = uint32(zeros(0,1));
            obj.Coordinates = zeros(0,obj.CoordinateSystem.Space.Dimension);
        end
        
        function newId = addPoints(obj, coordinates, ids, cSys)
            %addPoints add new points to the PointSet
            %   newId = addPoints(obj, coordinates, ids, cSys)
            
            % default coord sys is local
            if nargin < 4
                cSys = obj.CoordinateSystem;
            end
            
            % default id's are generated
            if nargin < 3 || isempty(ids)
                n = size(coordinates,1);
                ids = obj.generateIds(n);
            end
            
            % error check
            if numel(ids) ~= size(coordinates,1)
                error('NiHu:PointSet:invalid_argument', ...
                    'Number of Id arguments does not match number of coordinates');
            end
            if numel(unique(ids)) ~= numel(ids)
                error('NiHu:PointSet:invalid_argument', ...
                    'Id arguments are not unique');
            end
            if any(ismember(ids, obj.PointIds))
                error('NiHu:PointSet:invalid_argument', ...
                    'Id coincidence');
            end
            newId = ids(:);

            obj.PointIds = [obj.PointIds; newId];
            obj.Coordinates = [
                obj.Coordinates
                obj.CoordinateSystem.convertCoordinates(coordinates, cSys)
                ];
        end % of function addPoints
        
        function [coords, ids] = getCoordinatesByIndex(obj, idx, cSys)
            %getCoordinatesByIndex return coordinates by index
            %[coords, ids] = getCoordinatesByIndex(obj, idx, cSys)

            % default result coord sys is local
            if nargin < 3
                cSys = obj.CoordinateSystem;
            end
            
            % argument check
            if ~isa(cSys, 'CoordinateSystem')
                error('NiHu:PointSet:invalid_argument', ...
                    'argument ''cSys'' must be a ''CoordinateSystem'' instance');
            end
            
            % convert coordinates to result sys
            coords = cSys.convertCoordinates(obj.Coordinates(idx,:), obj.CoordinateSystem);
            
            % return ids if needed
            if nargout > 1
                ids = obj.PointIds(idx);
            end
        end
        
        function [coords, idx] = getCoordinatesById(obj, ids, cSys)
            %getCoordinatesByIndex return coordinates by index
            %[coords, ids] = getCoordinatesByIndex(obj, idx, cSys)

            % defautl csys is local
            if nargin < 3
                cSys = obj.CoordinateSystem;
            end
            
            % argument check
            if ~isa(cSys, 'CoordinateSystem')
                error('NiHu:PointSet:invalid_argument', ...
                    'argument ''cSys'' must be a ''CoordinateSystem'' instance');
            end
            
            idx = obj.id2idx(ids);
            coords = getCoordinatesByIndex(obj, idx, cSys);
        end % of function getCoordinatesById
        
        function newIds = merge(obj, otherSet)
            %merge merge other point set into this
            %newIds = merge(obj, otherSet)
            
            if ~isa(otherSet, 'PointSet')
                error('NiHu:PointSet:invalid_argument', ...
                    'Argument ''otherSet'' must be a ''PointSet'' instance');
            end
            
            newIds = otherSet.PointIds;
            
            % convert other coordinates to local system
            coords = obj.CoordinateSystem.convertCoordinates(otherSet.Coordinates, otherSet.CoordinateSystem);
            
            % check coincident IDs
            [ism, dst] = ismember(newIds, obj.PointIds);

            % add nodes with new IDs
            obj.addPoints(coords(~ism,:), otherSet.PointIds(~ism));

            % compare nodes with the same ID
            dst = dst(ism);
            src = find(ism);
            coinc = ismember(coords(src,:), obj.Coordinates(dst,:), 'rows');

            % add noncoincident nodes and ask for a new ID
            src = src(~coinc);
            newIds(src) = obj.addPoints(coords(src,:));
        end
    end
    
    methods (Access = private)
        function ids = generateIds(obj, n)
            if isempty(obj.PointIds)
                mx = uint32(0);
            else
                mx = max(obj.PointIds);
            end
            ids = mx + uint32(1 : n)';
        end
        
        function idx = id2idx(obj, id)
            [ism, idx] = ismember(id, obj.PointIds);
            if ~all(ism)
                error('NiHu:PointSet:invalid_argument', ...
                    'Id cannot found in PointSet');
            end
        end
    end
end
