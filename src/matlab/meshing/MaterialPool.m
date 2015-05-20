classdef MaterialPool < handle
    properties (SetAccess = private)
        ids
        materials
    end
    
    methods
        function obj = MaterialPool
            obj.ids = [];
            obj.materials = AcousticFluid.empty(1,0);
        end
        
        function id = add(obj, material, id)
            if ~isa(material, 'Material')
                error('NiHu:invalid_argument',...
                    'input should be a material');
            end
            
            if nargin < 3
                id = obj.assignNewId();
            else
                if ismember(obj.ids, id)
                    error('NiHu:out_of_range',...
                        'Unavailable material id: %d', id);
                end
            end
            obj.materials(end+1) = material;
            obj.ids(end+1) = id;
        end
    end % of public methods
    
    methods (Access = private)
        function id = assignNewId(obj)
            if ~isempty(obj.ids)
                id = max(obj.ids)+1;
            else
                id = 1;
            end
        end
    end % of private methods
end % of class
