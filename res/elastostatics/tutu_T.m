function tutu_T
nu = 1/3;
mu = 1;
r = [
    -.5 1 1 -1
    -1 -1 1 1
    0 0 0 0
    ];
c = mean(r,2);

N = size(r,2);

I = zeros(3,3);

for k = 1 : N
    tri = [c r(:,mod(k-1,N)+1), r(:,mod(k,N)+1)];
    T = rotmatrix(tri(:,1), tri(:,2), tri(:,3));
    tritri = T \ tri;
    tritri = bsxfun(@minus, tritri, tritri(:,1));
    R = tritri(1,2);
    theta2 = atan2(tritri(2,3), tritri(1,3));
    theta1 = atan2(tritri(2,2), tritri(1,2));
    I0 = primitive(mu, nu, R, theta2) - primitive(mu, nu, R, theta1);
    I = I + T*I0*T.';
end

disp(I)
end

function P = primitive(~, nu, R, theta)
a = sin(theta)*(1+log(R/cos(theta))) +...
    log( (cos(theta/2)-sin(theta/2))/(cos(theta/2)+sin(theta/2)) );
b = -cos(theta)*(1+log(R/cos(theta)));
P = R * [
    0, 0, a
    0, 0, b
    -a, -b, 0
    ] * (1-2*nu)/(8*pi*(1-nu));
end
