classdef Space
    %SPACE Definition of spaces
    
    properties (SetAccess = immutable)
        Id          % the space identified
        Dimension  % the space's dimensionality
    end
    
    methods (Access = private)
        function obj = Space(d)
            %SPACE Constructor
            % 
            obj.Id = uint32( d );
            obj.Dimension = d;
        end
    end
    
    enumeration
        D0(0)   % point space
        D1(1)   % 1D space
        D2(2)   % 2D space
        D3(3)   % 3D space
    end
end % of class
