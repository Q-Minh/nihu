N = 9;
[Xi, W] = duffy_line_line(N, 'edge');
f = @(x,y)log(x+y);
I = sum(f(Xi(:,1), Xi(:,2)) .* W);
syms x y;
I0 = double(int(int(f, x, 0, 1), y, 0, 1));
I/I0-1