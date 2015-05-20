classdef PropertyPool < handle
    properties (SetAccess = private)
        ids
        props
    end
    
    methods
        function obj = PropertyPool
            obj.ids = [];
            obj.props = AcousticVolume.empty(1,0);
        end
        
        function id = add(obj, prop, id)
            if ~isa(prop, 'Property')
                error('NiHu:invalid_argument',...
                    'input should be a property');
            end
            
            if nargin < 3
                id = obj.assignNewId();
            else
                if ismember(obj.ids, id)
                    error('NiHu:out_of_range',...
                        'Unavailable property id: %d', id);
                end
            end
            obj.props(end+1) = prop;
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
