function p = iemap(M,ps)
%IEMAP Coordinate mapping for infinite elements
if size(M,2) ~= size(ps,2);
    error('iemap3D:invalidArg','Mapping matrix has invalid size.');
end

switch (size(M,1))
    case 4 %Infinite quad element
        m = zeros(size(ps,1),4);
        s = ps(:,1);
        t = ps(:,2);
        m(:,1) = (1-t).*s./(s-1);
        m(:,2) = (1+t).*s./(s-1);
        m(:,3) = -(1-t).*(1+s)./(2*(s-1));
        m(:,4) = -(1+t).*(1+s)./(2*(s-1));
    case 6 %Infinite penta element
        m = zeros(size(ps,1),6);
        s = ps(:,1);
        t = ps(:,2);
        u = ps(:,3);
        m(:,1) =  (1-t).*s./(s-1);
        m(:,2) =  2*(1-(1+u)/2-(1-t)/2).*s./(s-1);
        m(:,3) =  (1+u).*s./(s-1);
        m(:,4) = -(1-t).*(1+s)./(2*(s-1));
        m(:,5) = -2*(1-(1+u)/2-(1-t)/2).*(1+s)./(2*(s-1));
        m(:,6) =  -(1+u).*(1+s)./(2*(s-1));
    case 8 %infinite hexa element
        m = zeros(size(ps,1),8);
        s = ps(:,1);
        t = ps(:,2);
        u = ps(:,3);
        m(:,1) =  1/2*(1-u).*(1-t).*s./(s-1);
        m(:,2) =  1/2*(1-u).*(1+t).*s./(s-1);
        m(:,3) =  1/2*(1+u).*(1+t).*s./(s-1);
        m(:,4) =  1/2*(1+u).*(1-t).*s./(s-1);
        m(:,5) = -1/2*(1-u).*(1-t).*(1+s)./(2*(s-1));
        m(:,6) = -1/2*(1-u).*(1+t).*(1+s)./(2*(s-1));
        m(:,7) = -1/2*(1+u).*(1+t).*(1+s)./(2*(s-1));
        m(:,8) = -1/2*(1+u).*(1-t).*(1+s)./(2*(s-1));
    otherwise
        error('iemap3D:invalidArg','Mapping matrix has invalid size.');
end
p = m*M;