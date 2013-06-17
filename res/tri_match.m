function [mu_lim, eta_lim] = tri_match(type)

switch type
    case 'edge'
        error('tria edge unimplemented');
    case 'face'
        eta_lim = [0 0; 1 0; 1 1];
        c = [1 0; 1 1; 0 1; -1 0; -1 -1; 0 -1];
        nc = size(c,1);
        for i = 1 : nc
            mu_lim{i} = [0 0; 0 0; c(i,:); c(mod(i,nc)+1,:)];
        end
    case 'corner'
        eta_lim = [-1 -1; 0 -1; 0 0];
        
        mu_lim{1} = [0 0; 0 0; -1 0; -1 -1];
        mu_lim{2} = [0 0; 0 0; -1 -1; 0 -1];

        mu_lim{3} = [-1 0; -1 0; -2 -1; -1 -1];
        mu_lim{4} = [0 -1; 0 -1; -1 -1; -1 -2];

        mu_lim{5} = [-1 -1; -2 -1; -2 -2; -1 -1];
        mu_lim{6} = [-1 -1; -2 -2; -1 -2; -1 -1];
end
