function [Xi, W] = duffy_general(mu_lim, eta_lim, N, inter_fun, map_fun)

Xi =  [];
W = [];

for q = 1 : length(mu_lim)
    [c_mu, w_mu] = quad_map(mu_lim{q}, N);
    
    for i = 1 : size(c_mu, 1)
        mu = c_mu(i,:);
        d = inter_fun(mu, eta_lim);
        [c_in, w_in] = map_fun(d, N);
        xi = [c_in [mu(1)+c_in(:,1), mu(2)+c_in(:,2)]];
        w = w_mu(i) * w_in;
        
        Xi = [Xi; xi];
        W = [W; w];
    end
end

end
