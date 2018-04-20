clear;

%%
q = 1/sqrt(3)/2;
c = ShapeSet.LinearQuad.Nodes / 2;
x0 = [-q -q];

%%

[Rvec, alphavec, thetavec, theta0vec] = plane_helper(c, x0);
nt = length(Rvec);

coeffs = [
    1 -1/(2*q) -1/(2*q) 1/(2*q)^2
    0 +1/(2*q) 0 -1/(2*q)^2
    0 0 0 1/(2*q)^2
    0 0 1/(2*q) -1/(2*q)^2
    ];

n = size(coeffs,1);
F = zeros(n,1);

for it = 1 : nt
    R = Rvec(it);
    theta = thetavec(it);
    a = alphavec(it);
    t0 = -theta0vec(it);
    
    for k = 1 : n
        
        a0 = coeffs(k,1);
        a1 = coeffs(k,2);
        a2 = coeffs(k,3);
        a12 = coeffs(k,4);
        
        b = R * sin(a);
        
        I = @(x) (a0 * cos(a + x)) / b ...
            + 2 * ( a1 * sin(a + t0) + a2 * cos(a + t0)) * atanh(cos(a) - sin(a) * tan(x/2)) ...
            - (a1 * sin(t0 - x) + a2 * cos(t0 - x)) * (log(b / sin(a + x)) + 1) ...
            + a12 / 2 * b * (sin(2 * (a + t0)) * log(cot((a + x)/2)) - 2 * sin(a + 2 * t0 - x));
        
        F(k) = F(k) + I(theta) - I(0);
        
    end
    
end

F / (4*pi)