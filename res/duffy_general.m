function [Xi, W] = duffy_general(N, type, match_fun, inter_fun, map_fun)

Xi =  [];
W = [];

[mu_lim, eta_lim] = match_fun(type);

for q = 1 : length(mu_lim)
    [c_mu, w_mu] = quad_map(mu_lim{q}, N);
    
    for i = 1 : size(c_mu, 1)
        mu = c_mu(i,:);
        d = inter_fun(mu, eta_lim);
        [c_in, w_in] = map_fun(d, N);
        xi = [c_in bsxfun(@plus, mu, c_in)];
        w = w_mu(i) * w_in;
        
        Xi = [Xi; xi];
        W = [W; w];
    end
end

end
