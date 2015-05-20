classdef Property < matlab.mixin.Heterogeneous & handle
    properties (Constant, GetAccess = protected)
        propClassIds = struct(...
            'AcousticVolume', 1,...
            'ElasticIsotropicSolid', 2);
    end % of constant protected properties
    
    properties (Abstract, Constant)
        classId     % unique class id
    end
    
    properties (SetAccess = immutable)
        name    % name of the property
    end
    
    methods
        function obj = Property(name)
            if ~ischar(name)
                error('NiHu:invalid_argument',...
                    'Property name should be a string');
            end
            obj.name = name;
        end
    end % of public methods
end
