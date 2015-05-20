classdef ElasticSolid < Material
    properties (Constant)
        classId = Material.matClassIds.ElasticSolid;
    end
    
    properties (SetAccess = immutable)
        massDensity     % mass density [kg/m^3]
        youngModulus    % Young's modulus of elasticity [Pa]
        poissonRatio    % Poisson's ratio [-]
    end
    
    properties (SetAccess = private)
        stiffnessMatrix     % 6x6 stiffness matrix
        complianceMatrix    % 6x6 compliance matrix
    end
    
    methods
        function obj = ElasticSolid(name, rho, E, nu)
            obj = obj@Material(name);
            
            if rho < 0
                error('NiHu:invalid_argument',...
                    'Mass density (%g) should not be negative', rho);
            end
            obj.massDensity = rho;
            
            obj.youngModulus = E;
            
            if nu < 0 || nu > .5
                error('NiHu:invalid_argument',...
                    'Posson''s ratio (%g) outside its range [0 0.5]', nu);
            end
            obj.poissonRatio = nu;
            
            obj.computeMatrices();
        end
    end
    
    methods (Access = private)
        function computeMatrices(obj)
            nu = obj.poissonRatio;
            G = obj.youngModulus/(2*(1+nu));
            lambda = obj.youngModulus*nu/(1+nu)/(1-2*nu);

            obj.stiffnessMatrix = [
                lambda*ones(3)+2*G*eye(3), zeros(3,3)
                zeros(3,3), G*eye(3)
                ];
            
            obj.complianceMatrix = 1/obj.youngModulus * [
                [1 -nu -nu; -nu 1 -nu; -nu -nu 1] zeros(3,3)
                zeros(3,3) 2*(1+nu)*eye(3)
                ];
        end
    end
    
end % of class
