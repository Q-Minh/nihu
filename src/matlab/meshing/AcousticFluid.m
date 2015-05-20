classdef AcousticFluid < Material
    %AcousticFluid Material describing an acoustic fluid characterised by
    %mass density and speed of sound.
    
    properties (Constant)
        classId = Material.matClassIds.AcousticFluid;
    end
    
    properties (SetAccess = immutable)
        massDensity     % mass density [kg/m^3]
        speedOfSound    % speed of sound [m/s]
    end
    
    methods
        function obj = AcousticFluid(name, rho, c)
            obj = obj@Material(name);
            if rho < 0
                error('NiHu:invalid_argument',...
                    'Mass density must not be negative');
            end
            obj.massDensity = rho;
            if c < 0
                error('NiHu:invalid_argument',...
                    'Speed of sound must not be negative');
            end
            obj.speedOfSound = c;
        end
    end % of public methods
    
end % of class AcousticFluid

