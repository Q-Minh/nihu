classdef Material < matlab.mixin.Heterogeneous & handle
    %MATERIAL Abstract base class for heterogeneous material collectors
    
    properties (Constant, GetAccess = protected)
        matClassIds = struct(...
            'AcousticFluid', 1,...
            'ElasticSolid', 2);
    end
    
    properties (Abstract, Constant)
        classId     % unique id of the class
    end
    
    properties (SetAccess = immutable)
        name    % user defined name of the material
    end
    
    methods
        function obj = Material(name)
            % only constructors can write names, we need to check here
            if ~ischar(name)
                error('NiHu:invalid_argument',...
                    'Material name should be a string');
            end
            obj.name = name;
        end
    end % of public methods
    
end % of class Material

