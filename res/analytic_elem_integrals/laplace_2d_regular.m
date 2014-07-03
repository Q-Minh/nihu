function I = laplace_2d_regular(type, x1, x2, y, a0, a1, a2)
switch type
    case 'slp'
        fun = @primitive_slp;
    case 'dlp'
        fun = @primitive_dlp;
    otherwise
        error('unknown kernel type: %s', type);
end
I = (fun(x1, y, a0, a1, a2) - fun(x2, y, a0, a1, a2))/(2*pi);
end

function I = primitive_slp(x, y, a0, a1, a2)
r2 = x^2 + y^2;
lnr = log(r2)/2;
I = (1/9 - lnr/3)*a2*x^3 + a1/2*(x^2/2 - lnr*r2) ...
    -a0*lnr*x + (a2*y^2/3-a0)*(y*atan(x/y)-x);
end

function I = primitive_dlp(x, y, a0, a1, a2)
r2 = x^2 + y^2;
lnr = log(r2)/2;
I = (a0-a2*y^2) * atan(x/y) + (a2*x+a1*lnr)*y;
end
