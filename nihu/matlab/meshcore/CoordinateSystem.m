classdef CoordinateSystem < handle
    properties (SetAccess = immutable)
        Space
    end
    
    properties (SetAccess = private)
        Origin
        Orientation
    end
    
    methods
        function obj = CoordinateSystem(space, origin, ori)
            if ~isa(space, 'Space')
                error('NiHu:CoordinateSystem:invalid_argument', ...
                    'Argument ''space'' needs to be a ''Space'' instance');
            end
            obj.Space = space;
            
            if numel(origin) ~= obj.Space.Dimension
                error('NiHu:CoordinateSystem:invalid_argument', ...
                    'Argument ''origin'' needs to be a %dD vector', obj.Space.Dimension);
            end
            obj.Origin = origin;
            
            if size(ori,2) ~= obj.Space.Dimension || size(ori,1) ~= obj.Space.Dimension-1
            end
            switch obj.Space
                case Space.D1 %#ok<PROP>
                    ori = 1;
                case Space.D2 %#ok<PROP>
                case Space.D3 %#ok<PROP>
                    ori(1,:) = ori(1,:) / norm(ori(1,:));
                    ori(2,:) = ori(2,:) - dot(ori(1,:), ori(2,:)) * ori(1,:);
                    ori(2,:) = ori(2,:) / norm(ori(2,:));
                    ori(3,:) = cross(ori(1,:), ori(2,:));
            end
            obj.Orientation = ori;
        end
        
        function translate(obj, v)
            if (~isnumeric(v) || any(size(v) ~= size(obj.Origin)))
                error('NiHu:CoordinateSystem:invalid_argument', ...
                    'Argument ''v'' must be a %dD vector', obj.Space.Dimension);
            end
            obj.Origin = obj.Origin + v;
        end
        
        function loc = convertCoordinates(obj, coords, cSys)
            if nargin < 3
                glob = bsxfun(@plus, coords * cSys.Orientation, cSys.Origin);
            else
                glob = coords;
            end
            loc = bsxfun(@minus, glob, obj.Origin) * obj.Orientation.';
        end
    end
end % of class CoordinateSystem
