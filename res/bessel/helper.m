function h = helper(nu, m)
if m == 0
    h = 1;
else
    h = helper(nu, m-1) * (4*nu^2-(2*m-1)^2)/(4*m);
end
