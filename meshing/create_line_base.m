function line = create_line_base(N)
% CREATE_LINE_BASE Creates basic line mesh. (NiHu / meshing)

x = linspace(-1, 1, N+1).';
    
line.Nodes = [(1:N+1).', x, zeros(size(x,1),2)];
line.Elements = [(1:N).', repmat([12 1 1], N, 1), (1:N).', (2:N+1).'];

% TODO: materials, properties?
line.Materials = [1 1 1 1 0 0];
line.Properties = [1 1 0 0 0 0];

end % create_line_base


