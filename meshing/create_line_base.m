function line = create_line_base(N)
%CREATE_LINE_BASE Create basic line mesh. (NiHu / meshing)

x = linspace(-1, 1, N+1).';
    
line.Nodes = [(1:N+1).', x, zeros(size(x,1),2)];
line.Elements = [(1:N).', repmat([12 1 1], N, 1), (1:N).', (2:N+1).'];

[line.Materials, line.Properties] = default_mat_prop();

end % create_line_base


