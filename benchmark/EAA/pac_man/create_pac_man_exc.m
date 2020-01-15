function [q_surf, p_in_field] = create_pac_man_exc(mesh, field, k, type)
%CREATE_PAC_MAN_EXC Create excitation for the Pac-Man problem
%
% Type:
%   'radiation' - Radiation from outer surface
%   'line_source' - Scattering from a line source
%
% References:
%   [1] https://eaa-bench.mec.tuwien.ac.at/fileadmin/t/eaa/PACMAN_benchmark_case.pdf

rho = 1.2041;
c = 343.21;

om = k * c;

N = size(mesh.Elements, 1);
q_surf = zeros(N,1);

switch type
    case {'rad', 'radiation'}
        if nargout > 1
            error('create_pac_man_excitation:invalid_output', ...
                'Only Q_SURF is available for radiation mode');
        end
        % Create constant field on the outer surface
        [~, sel] = mesh_select(mesh, 'abs(r-1) < 1e-4', 'ind');
        q_surf(sel) = 0.1 * (-1i*om*rho);
    case {'line', 'line_source'}
        % create line excitation
        r0 = [4 * cos(pi/4), 4 * sin(pi/4) 0];
        [x0, n0] = centnorm(mesh);
        [~, qs_in_line] = incident('line', r0, x0, n0, k);
        % The scattered field is the excitation
        q_surf = -qs_in_line * (4j);
        % Compute incident field in field points as needed
        if nargout > 1
            x0 = centnorm(field);
            p_in_field = incident('line', r0, x0, [], k);
        end
    otherwise
        error('create_pac_man_exc:unknown_type', ...
            'Unknown excitation type: %s', type);
end

end % of function create_pac_man_exc
