function D = ibemD(model, k, type, points, pairs)

%% parameter check
if nargin < 5
    pairs = [];
end
if nargin < 4
    points = [];
end

%% Gaussian quadrature divisions
% singular near mid far
gauss3 = assemble_gauss_struct(3, [9 7 5 2]);
gauss4 = assemble_gauss_struct(4, [9 7 5 2]);
[nodes, elements] = extract_bem_mesh(model);
dist = nodes(elements(:,3)+1,:) - nodes(elements(:,2)+1,:);
dist = [2 5] * max(sqrt(dot(dist,dist,2)));

%% System matrices
if ~isempty(pairs)
    error('The sparse option of ibemB is not yet implemented.');
else
    %% full matrices (BEM mode)
    if ~isempty(points) % full field point matrices
        switch lower(type)
            case 'const'
                [D, ~] = bemHG_const(nodes, elements, gauss3, gauss4, dist, k, points);
            case 'lin'
                [D, ~] = bemHG_lin(nodes, elements, gauss3, gauss4, dist, k, points);
        end
    else % full surface matrices
        switch lower(type)
            case 'const'
                D = ibemD_const(nodes, elements, gauss3, gauss4, dist, k);
            case 'lin'
                D = ibemD_lin(nodes, elements, gauss3, gauss4, dist, k);
        end
    end
end
end