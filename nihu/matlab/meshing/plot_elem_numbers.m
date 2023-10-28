function plot_elem_numbers(model)
%PLOT_ELEM_NUMBERS Plot element normals of a mesh
%   PLOT_ELEM_NORMALS(MESH)

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

% Last modified 02.12.2009.

%% Argument check
narginchk(1, 1);

%% plot lines and markers
[start, gnorm, w] = geo2gauss(model, 1);
stop = start + bsxfun(@times, sqrt(w), gnorm);
% line([start(:,1) stop(:,1)].',...
%     [start(:,2) stop(:,2)].',...
%     [start(:,3) stop(:,3)].',...
%     'Color', 'black');
% hold on;
% plot3(start(:,1), start(:,2), start(:,3), '.k');
text(stop(:,1), stop(:,2), stop(:,3), num2str(model.Elements(:,1)));
% hold off;
end
