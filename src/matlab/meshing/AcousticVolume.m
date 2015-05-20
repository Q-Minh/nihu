classdef AcousticVolume < Property
    properties (Constant)
        classId = Property.propClassIds.AcousticVolume;
    end
    
    properties (SetAccess = immutable)
        material    % the acoustic material
    end
    
    methods
        function obj = AcousticVolume(name, material)
            obj = obj@Property(name);
            if ~isa(material, 'AcousticFluid')
                error('NiHu:invalid_argument',...
                    'Material should be an acoustic fluid');
            end
            obj.material = material;
        end
    end % of public methods
end % of class
