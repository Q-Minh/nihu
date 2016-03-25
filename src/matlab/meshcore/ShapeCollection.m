classdef ShapeCollection < handle
    methods (Static = true)
        function obj = get_instance(id, varargin)
            persistent ShapeFunctions
            if isempty(ShapeFunctions)
                ShapeFunctions = cell(1,1);
            end
            if numel(ShapeFunctions) < id || isempty(ShapeFunctions{id})
                ShapeFunctions{id} = ShapeSet(id, varargin{:});
            end
            obj = ShapeFunctions{id};
        end % of function
    end % of static methods
end % of class
