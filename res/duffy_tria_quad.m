function [Xi, W] = duffy_tria_quad(N, type)

switch lower(type)
    case 'edge'
        [xi1, w1] = duffy_tria_tria(N, 'edge');
        [Xi1, W1] = quadrature_map(xi1(:,3:4), w1, [-1 -1; 1 -1; 1 1], 23);
        [xi2, w2] = duffy_tria_tria(N, 'corner');
        [Xi2, W2] = quadrature_map(xi2(:,3:4), w2, [-1 -1; 1 1; -1 1], 23);
        Xi = [xi1(:,1:2) Xi1; xi2(:,1:2) Xi2];
        W = [W1; W2];
    case 'corner'
        [xi, w] = duffy_tria_tria(N, 'corner');
        [Xi1, W1] = quadrature_map(xi(:,3:4), w, [-1 -1; 1 -1; 1 1], 23);
        [Xi2, W2] = quadrature_map(xi(:,3:4), w, [-1 -1; 1 1; -1 1], 23);
        Xi = [xi(:,1:2) Xi1; xi(:,1:2) Xi2];
        W = [W1; W2];
end

end
