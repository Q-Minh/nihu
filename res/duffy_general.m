function [Xi, W] = duffy_general(mu_lim, Teta, N)

eta_lim = [
    0 0;
    1 0;
    1 1;
    0 1
    ];

xi_lim = [
    0 0;
    1 0;
    1 1;
    0 1
    ];

Xi =  [];
W = [];

for q = 1 : length(mu_lim)
    % erre duffintani
    [c_mu, w_mu] = quad_map(mu_lim{q}, N);
    
    for i = 1 : size(c_mu, 1)
        mu = c_mu(i,:);
        eta = bsxfun(@plus, mu, -eta_lim*Teta);
        d = quad_intersect(eta, xi_lim);
        [c_in, w_in] = quad_map(d, N);
        xi = [c_in [mu(1)-c_in(:,1), mu(2)-c_in(:,2)]/Teta];
        w = w_mu(i) * w_in;
        
        Xi = [Xi; xi];
        W = [W; w];
    end
end

end
