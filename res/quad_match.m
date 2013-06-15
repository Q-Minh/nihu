function [mu_lim, eta_lim] = quad_match(type)

switch type
    case 'edge'
        eta_lim = [0 -1; 1 -1; 1 0; 0 0];

        mu_lim{1} = [0 0; 0 0; 1 0; 1 -1];
        mu_lim{2} = [0 0; 0 0; 1 -1; 0 -1];
        
        mu_lim{3} = [
            0 0;
            0 0;
            0 -1;
            -1 -1
            ];
        
        mu_lim{4} = [
            0 0
            0 0
            -1 -1
            -1 0
            ];
        
        mu_lim{5} = [
            -1 -1
            0 -1
            0 -2
            -1 -2
            ];
        
        mu_lim{6} = [
            0 -1
            1 -1
            1 -2
            0 -2
            ];
    case 'face'
        c = [
            1 0
            1 1
            0 1
            -1 1
            -1 0
            -1 -1
            0 -1
            1 -1
            ];
        for i = 1 : 8
            mu_lim{i} = [0 0; 0 0; c(i,:); c(mod(i,8)+1,:)];
        end
        eta_lim = [0 0; 1 0; 1 1; 0 1];
    case 'corner'
        mu_lim{1} = [0 0; 0 0; -1 0; -1 -1];
        mu_lim{2} = [0 0; 0 0; -1 -1; 0 -1];

        mu_lim{3} = [-1 0; -2 0; -2 -1; -1 -1];
        mu_lim{4} = [-1 -1; -2 -1; -2 -2; -1 -2];
        mu_lim{5} = [0 -1; 0 -2; -1 -2; -1 -1];
        
        eta_lim = [-1 -1; 0 -1; 0 0; -1 0];
end
