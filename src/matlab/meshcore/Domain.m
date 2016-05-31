classdef Domain
    %DOMAIN Intrinsic domain representations
    
    % Last modified: 2015.03.31. FP. Point domain introduced
    
    properties (SetAccess = immutable)
        Id
        Space
        CornerNodes
        Volume
        Center
    end
    
    methods (Access = private)
        function obj = Domain(space, cornernodes, volume, center)
            % DOMAIN constructor
            % obj = Domain(space, cornernodes, volume, center)
            
            if ~isa(space, 'Space')
                error('NiHu:Domain:invalid_argument', ...
                    'Argument ''space'' must be a ''Space'' instance');
            end
            
            if (size(cornernodes,2) ~= space.Dimension)
                error('NiHu:Domain:invalid_argument', ...
                    'Corner node - space dimension mismatch');
            end
            
            obj.Id = uint32( space.Id * 10 + size(cornernodes,1) );
            
            obj.Space = space;
            obj.CornerNodes = cornernodes;
            obj.Volume = volume;
            obj.Center = center;
        end
    end % of private methods
    
    enumeration
        Point(Space.D0, zeros(1,0), 0, zeros(1,0));
        Line(Space.D1, [-1; +1], 2, 0);
        Tria(Space.D2, [0 0; 1 0; 0 1], .5, [1 1]/3);
        Quad(Space.D2, [-1 -1; 1 -1; 1 1; -1 1], 4, [0 0]);
        Tetra(Space.D3, [0 0 0; 1 0 0; 0 1 0; 0 0 1], 1/6, [1 1 1]/4);
        Penta(Space.D3, [0 0 -1; 1 0 -1; 0 1 -1; 0 0 1; 1 0 1; 0 1 1], 1, [1 1 0]/3);
        Hexa(Space.D3, [
            -1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1
            -1 -1 +1; 1 -1 +1; 1 1 +1; -1 1 +1
            ], 8, [0 0 0]);
    end
end
