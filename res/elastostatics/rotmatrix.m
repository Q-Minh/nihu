function T = rotmatrix(r0, r1, r2)
e3 = cross(r1-r0, r2-r0);
e3 = e3/norm(e3);
e2 = r2-r1;
e2 = e2/norm(e2);
e1 = cross(e2, e3);
T = [e1 e2 e3];
end