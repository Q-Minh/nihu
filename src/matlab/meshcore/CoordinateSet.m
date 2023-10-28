classdef CoordinateSet < handle
    %COORDINATESET A set of coordinates
    
    properties (SetAccess = private)
        coordinates % coordinates of a point in a row
    end
    
    methods
        function obj = CoordinateSet(varargin)
            switch nargin
                case 0
                case 1
                    xi = varargin{1};
                    obj.coordinates = xi;
            end
        end
        
        function nC = getNumCoordinates(obj)
            nC = size(obj.coordinates,1);
        end
        
        function nD = getDimension(obj)
            nD = size(obj.coordinates,2);
        end
        
        function [ism, idx] = find(obj, coords, tol)
            [ism, idx] = ismember(round(coords/tol), round(obj.coordinates/tol), 'rows');
        end
        
        function idx = addCoordinates(obj, coords, method, tol)
            switch method
                case 'merge'
                    [ism, where] = obj.find(coords, tol);
                    idx(ism) = where(ism);
                    idx(~ism) = obj.getNumCoordinates + (1:sum(~ism));
                    idx = uint32(idx(:));
                    obj.coordinates = [
                        obj.coordinates
                        coords(~ism, :)
                        ];
                case 'concat'
                    idx = uint32((obj.getNumCoordinates + (1:size(coords,1))).');
                    obj.coordinates(idx,:) = coords;
                otherwise
                    % error here
            end
        end
        
        function [ia, ic] = mergeCoincident(obj, tol)
            [~, ia, ic] = unique(round(obj.coordinates/tol), 'rows', 'stable');
            obj.coordinates = obj.coordinates(ia,:);
        end
    end
end

