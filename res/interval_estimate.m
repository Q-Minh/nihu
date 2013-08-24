clear;
clc;

y = sym('y', 'real');       % distance on the element
D = sym('D', 'positive');   % elemen size
k = sym('k', 'positive');   % relative distance from source

r = k*D + y;    % absolute distance from source

n_max = 3;      % maximal inverse power r^(-n)
kvec = (1:.1:20)';  % relative distance scale

% preallocate for result
order = zeros(length(kvec), n_max);


N = 15;     % maximal polynomial order
eps = 1e-2; % prescribed error

for n = 1 : n_max % for each green power
    disp(n);
    
    g = 1/r^n;  % Green's function

    % Taylor expansion on the element
    Taylor = sym('T', [1, N+1]);
    for i = 0 : N
        Taylor(i+1) = diff(g, y, i)/factorial(i);
    end
    Taylor = simple(subs(Taylor, y, 0));
    Taylor = simple(Taylor .* y.^(0:N));
    % integrate the taylor expansion over the element
    I = simple(int(Taylor, y, -D/2, D/2));
    % integrate the Green's function over the element
    Ianal = simple(int(g, y, -D/2, D/2));
    
    % compute relative error for each order and each distance
    cont = sym('cont', [length(kvec), N+1]);
    for kk = 1 : length(kvec)
        cont(kk,:) = subs(I/Ianal, k, kvec(kk));
    end
    % cumulated relative error
    Error = cumsum(double(cont),2) - 1;
    % compute minimal order
    order(:,n) = sum(abs(Error) > eps, 2);
end

plot(kvec, order);
xlabel('relative distance');
ylabel('order');
