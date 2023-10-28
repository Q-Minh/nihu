classdef PointSet < handle
    %PointSet
    %   class representing a set of coordinates
    
    properties (SetAccess = private)
        CoordinateSystem
        PointIds
        Coordinates
    end % of private properties
    
    methods
        function obj = PointSet(varargin)
            %PointSet constructor
            % obj = PointSet(id, coordSys)
            
            if (isnumeric(varargin{1}))
                % 
                
            
            elseif isa(varargin{1}, 'CoordinateSystem')
                % 
                coordSys = varargin{1};
            else
                error('NiHu:PointSet:invalid_argument', ...
                    'First argument must be an array of coordinates or ''CoordinateSystem'' instance');
            end
            
            obj.CoordinateSystem = coordSys;
            
            % allocate PointIds and Coordinates
            obj.PointIds = uint32(zeros(0,1));
            obj.Coordinates = zeros(0,obj.CoordinateSystem.Space.Dimension);
        end
        
        function newId = addPoints(obj, coordinates, ids, cSys)
            %addPoints add new points to the PointSet
            %   newId = addPoints(obj, coordinates, ids, cSys)
            %    the points are transformed from the source to the local
            %    coordinate system. Id coincidence causes an error.
            
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
        end % of function getCoordinatesByIndex
        
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
        
        function newIds = merge(obj, otherSet, tol)
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
            dst = dst(ism);
            src = find(ism);

            % add nodes with new IDs
            obj.addPoints(coords(~ism,:), otherSet.PointIds(~ism));
            
            % If tolerance parameter was given
            if nargin > 2 
            
                % skip close point pairs
                diff = coords(src,:) - obj.Coordinates(dst,:);
                skip = dot(diff,diff,2) < tol^2;

                % add far points and ask for new id
                newIds(src(~skip)) = obj.addPoints(coords(src(~skip),:));
            else
                newIds(src) = obj.addPoints(coords(src, :));
            end
        end % of function merge
        
        function translate(obj, v)
            % TRANSLATE Translate the coordinate system of the PointSet
            obj.CoordinateSystem.translate(v);
        end
        
        
    end % of methods
    
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
                    'Some Ids not found in PointSet');
            end
        end
    end
end % of class PointSet
