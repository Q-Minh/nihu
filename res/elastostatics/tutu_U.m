function tutu_U
nu = 1/3;
mu = 1;
r = [
    -.5 1 1 -1
    -1 -1 1 1
    0 0 0 0
    ];
c = mean(r,2);

n = size(r,2);

I = zeros(3,3);

for k = 1 : n
    tri = [c r(:,mod(k-1,n)+1), r(:,mod(k,n)+1)];
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

function P = primitive(mu, nu, R, theta)
q = 2*atanh(tan(theta/2));
P = R * [
    sin(theta)+(3-4*nu)*q, -(cos(theta)+1), 0
    -(cos(theta)+1), 4*(1-nu)*q-sin(theta), 0
    0, 0, (3-4*nu)*q
    ] / (16*pi*mu*(1-nu));
end
